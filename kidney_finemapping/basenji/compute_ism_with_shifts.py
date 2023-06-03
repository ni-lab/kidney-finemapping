#!/usr/bin/env python
import json
import os
from optparse import OptionParser

import h5py
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from kidney_finemapping.basenji.basenji_utils import seqnn
from kidney_finemapping.basenji.basenji_utils import stream
from kidney_finemapping.basenji.basenji_utils import vcf as bvcf


def main():
    """
    Run in silico saturation mutagenesis on a set of SNPs where the ISM scores are determined by taking the mean
    across the 1344 shifts in the receptive field.
    - args[0] <params_file> - JSON file specifying model parameters
    - args[1] <model_file> - Model file
    - args[2] <vcf_file> - VCF file
    """
    usage = "usage: %prog [options] <params_file> <model_file> <vcf_file>"
    parser = OptionParser(usage)
    parser.add_option("-f", dest="genome_fasta",
                      default=None,
                      help="Genome FASTA for sequences [Default: %default]")
    parser.add_option("-o", dest="out_dir",
                      default="sad",
                      help="Output directory for tables and plots [Default: %default]")
    parser.add_option("-l", dest="mut_len",
                      default=20, type="int",
                      help="Length of center sequence to mutate. Default is 10 upstream and 10 downstream of variant. [Default: %default]")
    parser.add_option("--rc", dest="rc",
                      default=False, action="store_true",
                      help="Average forward and reverse complement predictions [Default: %default]")
    parser.add_option("--shifts", dest="shifts",
                      default="0", type="str",
                      help="Ensemble prediction shifts [Default: %default]")
    parser.add_option("-t", dest="targets_file",
                      default=None, type="str",
                      help="File specifying target indexes and labels in table format")
    (options, args) = parser.parse_args()
    # Setup
    num_expected_args = 3
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    # Parse arguments
    params_file = args[0]
    model_file = args[1]
    vcf_file = args[2]

    os.makedirs(options.out_dir, exist_ok=True)
    options.shifts = [int(shift) for shift in options.shifts.split(",")]

    assert (options.mut_len > 0)
    options.mut_up = options.mut_len // 2
    options.mut_down = options.mut_len - options.mut_up

    # Read model parameters and model setup
    with open(params_file) as params_open:
        params = json.load(params_open)
    params_model = params["model"]
    params_train = params["train"]

    targets_df = pd.read_csv(options.targets_file, sep="\t", index_col=0)
    target_slice = targets_df.index

    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(model_file)
    seqnn_model.build_slice(target_slice)
    seqnn_model.build_ensemble(options.rc, options.shifts)

    num_targets = seqnn_model.num_targets()

    #################################################################
    # load SNPs
    snps = bvcf.vcf_snps(vcf_file, validate_ref_fasta=options.genome_fasta, flip_ref=True)  # snps are flipped here so that ref matches reference genome

    # open genome FASTA
    genome_open = pysam.Fastafile(options.genome_fasta)

    def seq_gen(full_seq, mut_i, validation_ism_seq):
        center_full_seq = params_model["seq_length"] + mut_i - 1
        assert full_seq[center_full_seq] == validation_ism_seq[mut_i], "Center of full seq should match corresponding position in ISM seq"

        pos_shifts_seq = full_seq[center_full_seq - (params_model["seq_length"] - 1):
                                  center_full_seq + params_model["seq_length"]]

        # Create ISM mutation
        for allele in ["A", "C", "G", "T"]:
            center = len(pos_shifts_seq) // 2
            assert pos_shifts_seq[center] == validation_ism_seq[mut_i], "Center of pos shifts seq should match correpsonding position in ISM seq"
            mut_pos_shifts_seq = pos_shifts_seq[:center] + allele + pos_shifts_seq[center + 1:]

            # For each ISM mutation, iterate over all pos shifts
            for i in range(params_model["seq_length"]):
                seq = mut_pos_shifts_seq[i:i + params_model["seq_length"]]
                seq_1hot, _ = bvcf.dna_length_1hot(seq, params_model["seq_length"])
                yield seq_1hot

    ref_ism_preds_dir = f"{options.out_dir}/ref_ism_preds"
    alt_ism_preds_dir = f"{options.out_dir}/alt_ism_preds"

    os.makedirs(ref_ism_preds_dir, exist_ok=True)
    os.makedirs(alt_ism_preds_dir, exist_ok=True)

    ref_ism_seqs = []
    alt_ism_seqs = []
    snp_ids = []
    snp_flips = [snp.flipped for snp in snps]

    for snp in tqdm(snps):
        snp_ids.append(snp.rsid)

        ref_ism_preds_file = f"{ref_ism_preds_dir}/{snp.rsid}"
        alt_ism_preds_file = f"{alt_ism_preds_dir}/{snp.rsid}"

        # REFERENCE
        # Assume symmetric mutation range to avoid dealing with strandedness
        mut_start = snp.pos - options.mut_up  # 1-based
        mut_end = snp.pos + len(snp.ref_allele) + options.mut_down  # 1-based, exclusive

        # Get ISM region
        ref_ism_seq = genome_open.fetch(snp.chr, mut_start - 1, mut_end - 1)  # 0-based, exclusive endpoint
        assert len(ref_ism_seq) == len(snp.ref_allele) + options.mut_len, "ISM region does not have the correct length"
        assert ref_ism_seq[options.mut_up:options.mut_up + len(snp.ref_allele)] == snp.ref_allele, "Genome does not match SNP ref allele"
        ref_ism_seqs.append(ref_ism_seq)

        # Get full sequence range for pos shifts
        full_seq_start = mut_start - (params_model["seq_length"] - 1)  # 1-based
        full_seq_end = mut_end + (params_model["seq_length"] - 1)  # 1-based, exclusive
        ref_full_seq = genome_open.fetch(snp.chr, full_seq_start - 1, full_seq_end - 1)  # 0-based, exclusive endpoint
        assert len(ref_full_seq) == len(snp.ref_allele) + options.mut_len + 2 * (params_model["seq_length"] - 1), \
                "Full sequence does not have the correct length"

        ref_ism_preds = []
        for mut_i in range(len(ref_ism_seq)):
            mut_i_pos_shift_preds = np.zeros((4, params_model["seq_length"], num_targets))
            preds_stream = stream.PredStreamGen(seqnn_model,
                                                seq_gen(ref_full_seq, mut_i, validation_ism_seq=ref_ism_seq),
                                                params_train["batch_size"])

            for pi in range(4 * params_model["seq_length"]):
                mut_i_pos_shift_preds[pi // params_model["seq_length"], pi % params_model["seq_length"]] = preds_stream[pi]

            mut_i_preds = mut_i_pos_shift_preds.mean(axis=1)  # take average across pos shifts
            ref_ism_preds.append(mut_i_preds)

        ref_ism_preds = np.array(ref_ism_preds)
        np.save(ref_ism_preds_file, ref_ism_preds)

        # ALTERNATE
        assert len(snp.alt_alleles) == 1, "Each SNP should only have one alt allele"
        alt_allele = snp.alt_alleles[0]

        # Assume symmetric mutation range to avoid dealing with strandedness
        mut_start = snp.pos - options.mut_up  # 1-based
        mut_end = snp.pos + len(snp.ref_allele) + options.mut_down  # 1-based, exclusive

        # Get ISM region
        alt_ism_seq = genome_open.fetch(snp.chr, mut_start - 1, mut_end - 1)  # 0-based, exclusive endpoint
        alt_ism_seq = alt_ism_seq[:options.mut_up] + alt_allele + alt_ism_seq[options.mut_up + len(snp.ref_allele):]  # change to alternate sequence
        assert len(alt_ism_seq) == len(alt_allele) + options.mut_len, "ISM region does not have the correct length"
        assert alt_ism_seq[options.mut_up:options.mut_up + len(alt_allele)] == alt_allele, "Alt allele doesn't match ism seq"
        alt_ism_seqs.append(alt_ism_seq)

        # Get full sequence range for pos shifts
        full_seq_start = mut_start - (params_model["seq_length"] - 1)  # 1-based
        full_seq_end = mut_end + (params_model["seq_length"] - 1)  # 1-based, exclusive
        alt_full_seq = genome_open.fetch(snp.chr, full_seq_start - 1, full_seq_end - 1)  # 0-based, exclusive endpoint
        alt_full_seq = alt_full_seq[:params_model["seq_length"] - 1 + options.mut_up] + alt_allele + \
                       alt_full_seq[params_model["seq_length"] - 1 + options.mut_up + len(snp.ref_allele):]  # change to alt sequence

        assert len(alt_full_seq) == len(alt_allele) + options.mut_len + 2 * (params_model["seq_length"] - 1), \
                "Full sequence does not have the correct length"
        assert alt_full_seq[params_model["seq_length"] - 1:params_model["seq_length"] - 1 + len(alt_ism_seq)] == alt_ism_seq, "alt_full_seq center does not match alt_ism_seq"

        alt_ism_preds = []
        for mut_i in range(len(alt_ism_seq)):
            mut_i_pos_shift_preds = np.zeros((4, params_model["seq_length"], num_targets))
            preds_stream = stream.PredStreamGen(seqnn_model,
                                                seq_gen(alt_full_seq, mut_i, validation_ism_seq=alt_ism_seq),
                                                params_train["batch_size"])

            for pi in range(4 * params_model["seq_length"]):
                mut_i_pos_shift_preds[pi // params_model["seq_length"], pi % params_model["seq_length"]] = preds_stream[pi]

            mut_i_preds = mut_i_pos_shift_preds.mean(axis=1)  # take average across pos shifts
            alt_ism_preds.append(mut_i_preds)

        alt_ism_preds = np.array(alt_ism_preds)
        np.save(alt_ism_preds_file, alt_ism_preds)

    # close genome
    genome_open.close()

    # Write rsids, ref_ism_seqs, and alt_ism_seqs to h5
    info_h5 = f"{options.out_dir}/info.h5"
    with h5py.File(info_h5, "w") as h5_out:
        h5_out.create_dataset("rsid", data=np.array(snp_ids, dtype="S"))
        h5_out.create_dataset("ref_ism_seq", data=np.array(ref_ism_seqs, dtype="S"))
        h5_out.create_dataset("alt_ism_seq", data=np.array(alt_ism_seqs, dtype="S"))
        h5_out.create_dataset("snp_flipped", data=snp_flips)


if __name__ == "__main__":
    main()
