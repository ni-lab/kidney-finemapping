#!/usr/bin/env python
from __future__ import print_function

import json
import os
import pickle
from optparse import OptionParser
from queue import Queue
from threading import Thread

import h5py
import numpy as np
import pandas as pd
import pysam

from kidney_finemapping.basenji.basenji_utils import seqnn
from kidney_finemapping.basenji.basenji_utils import stream
from kidney_finemapping.basenji.basenji_utils import vcf as bvcf


def main():
    usage = 'usage: %prog [options] <params_file> <model_file> <vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='genome_fasta',
                      default='%s/data/hg19.fa' % os.environ['BASENJIDIR'],
                      help='Genome FASTA for sequences [Default: %default]')
    parser.add_option('--rc', dest='rc',
                      default=False, action='store_true',
                      help='Average forward and reverse complement predictions [Default: %default]')
    parser.add_option('--shifts', dest='shifts',
                      default='0', type='str',
                      help='Ensemble prediction shifts [Default: %default]')
    parser.add_option('-t', dest='targets_file',
                      default=None, type='str',
                      help='File specifying target indexes and labels in table format')
    parser.add_option('-o', dest='out_dir',
                      default='sad',
                      help='Output directory for tables and plots [Default: %default]')
    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 3
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    params_file = args[0]
    model_file = args[1]
    vcf_file = args[2]

    os.makedirs(options.out_dir, exist_ok=True)

    options.shifts = [int(shift) for shift in options.shifts.split(',')]

    # read parameters and targets
    # read model parameters
    with open(params_file) as params_open:
        params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

    targets_df = pd.read_csv(options.targets_file, sep='\t', index_col=0)
    target_ids = targets_df.identifier
    target_labels = targets_df.description
    target_slice = targets_df.index

    # setup model
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(model_file)
    seqnn_model.build_slice(target_slice)
    seqnn_model.build_ensemble(options.rc, options.shifts)

    targets_length = seqnn_model.target_lengths[0]

    # load SNPs
    # read SNPs form VCF
    snps = bvcf.vcf_snps(vcf_file)  # does not flip SNPs according to ref fasta

    num_snps = len(snps)
    num_pos = params_model['seq_length']

    # open genome FASTA
    genome_open = pysam.Fastafile(options.genome_fasta)

    def snp_gen():
        for snp in snps:
            for i in range(num_pos):
                # get SNP sequences
                snp_1hot_list = bvcf.snp_seq1_pos_shifts(snp, params_model['seq_length'], i, genome_open)
                for snp_1hot in snp_1hot_list:
                    yield snp_1hot

    #################################################################
    # setup output

    sad_out = initialize_output_h5(options.out_dir, snps, target_ids, target_labels, targets_length)

    # predict SNP scores, write output
    # initialize predictions stream
    preds_stream = stream.PredStreamGen(seqnn_model, snp_gen(), params_train['batch_size'])

    # predictions index
    pi = 0

    for si in range(num_snps * num_pos):
        # get predictions
        ref_preds = preds_stream[pi]
        pi += 1
        alt_preds = preds_stream[pi]
        pi += 1

        # process SNP
        write_snp(ref_preds, alt_preds, sad_out, si)

    # close files
    genome_open.close()
    sad_out.close()


def initialize_output_h5(out_dir, snps, target_ids, target_labels, targets_length):
    """Initialize an output HDF5 file for SAD stats."""

    num_targets = len(target_ids)
    num_snps = len(snps)
    num_pos = 1344

    sad_out = h5py.File('%s/sad.h5' % out_dir, 'w')

    # write SNPs
    snp_ids = np.array([snp.rsid for snp in snps], 'S')
    sad_out.create_dataset('snp', data=snp_ids)

    # write SNP chr
    snp_chr = np.array([snp.chr for snp in snps], 'S')
    sad_out.create_dataset('chr', data=snp_chr)

    # write SNP pos
    snp_pos = np.array([snp.pos for snp in snps], dtype='uint32')
    sad_out.create_dataset('pos', data=snp_pos)

    # check flips
    snp_flips = [snp.flipped for snp in snps]

    # write SNP reference allele
    snp_refs = []
    snp_alts = []
    for snp in snps:
        if snp.flipped:
            snp_refs.append(snp.alt_alleles[0])
            snp_alts.append(snp.ref_allele)
        else:
            snp_refs.append(snp.ref_allele)
            snp_alts.append(snp.alt_alleles[0])
    snp_refs = np.array(snp_refs, 'S')
    snp_alts = np.array(snp_alts, 'S')
    sad_out.create_dataset('ref_allele', data=snp_refs)
    sad_out.create_dataset('alt_allele', data=snp_alts)

    # write targets
    sad_out.create_dataset('target_ids', data=np.array(target_ids, 'S'))
    sad_out.create_dataset('target_labels', data=np.array(target_labels, 'S'))

    # initialize SAD stats
    for sad_stat in ["REF", "ALT"]:
        sad_out.create_dataset(sad_stat,
                               shape=(num_snps * num_pos, targets_length, num_targets),
                               dtype='float16')

    sad_out.create_dataset("SAD",
                           shape=(num_snps * num_pos, num_targets),
                           dtype='float16')

    return sad_out


def write_snp(ref_preds, alt_preds, sad_out, si):
    """Write SNP predictions to HDF."""
    ref_preds = ref_preds.astype('float64')
    alt_preds = alt_preds.astype('float64')

    # sum across length
    ref_preds_sum = ref_preds.sum(axis=0)
    alt_preds_sum = alt_preds.sum(axis=0)

    # compare reference to alternative via mean subtraction
    sad = alt_preds_sum - ref_preds_sum
    sad_out['SAD'][si, :] = sad.astype('float16')

    # predictions
    sad_out['REF'][si, :] = ref_preds.astype('float16')
    sad_out['ALT'][si, :] = alt_preds.astype('float16')


if __name__ == '__main__':
    main()
