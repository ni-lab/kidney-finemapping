#!/usr/bin/env python
import glob
import os
from optparse import OptionParser

import h5py
import numpy as np
import pandas as pd
from natsort import natsorted


def main():
    """
    Merge SAD h5 files from output of compute_sad.py from each chromosome into a single h5 file.
    - args[0] = sad_chr_dir: Directory containing SAD directories from output of compute_sad.py for each chromosome.
    """
    usage = "usage: %prog [options] <sad_chr_dir>"
    parser = OptionParser(usage)
    parser.add_option("-n", dest="num_files",
                      default=22, type=int,
                      help="Expected number of files to merge. [Default: %default]")
    parser.add_option("--vcf", dest="vcf",
                      default=False, action="store_true",
                      help="If set, combine SNP VCFs from SAD output to out directory as well.")
    parser.add_option("-o", dest="out_dir", default="all_chrs", help="Output directory for merged SAD h5 file.")
    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    os.makedirs(options.out_dir, exist_ok=True)

    # Parse args
    sad_dir = args[0]

    # Create h5 containing SADs from all batches
    sad_files = natsorted(glob.glob(os.path.join(sad_dir, "chr*/sad.h5")))
    assert len(sad_files) == options.num_files, f"Expected {options.num_files} sad files for {options.num_files} " \
                                                f"chromosomes, got {len(sad_files)} files"

    sads = [h5py.File(sad_file, mode="r") for sad_file in sad_files]
    sad_h5_file = os.path.join(options.out_dir, "sad.h5")
    sad_h5 = h5py.File(sad_h5_file, "w")

    sad_h5.create_dataset("target_ids", data=np.array(sads[0]["target_ids"]))
    sad_h5.create_dataset("target_labels", data=np.array(sads[0]["target_labels"]))
    sad_h5.create_dataset("snp", data=np.array([snp for sad in sads for snp in sad["snp"]], dtype="S"))
    sad_h5.create_dataset("ref_allele", data=np.array([allele for sad in sads for allele in sad["ref_allele"]], dtype="S"))
    sad_h5.create_dataset("alt_allele", data=np.array([allele for sad in sads for allele in sad["alt_allele"]], dtype="S"))
    sad_h5.create_dataset("snp_flipped", data=np.array([pos for sad in sads for pos in sad["snp_flipped"]]))
    sad_h5.create_dataset("chr", data=np.array([chr_num for sad in sads for chr_num in sad["chr"]], dtype="S"))
    sad_h5.create_dataset("pos", data=np.array([pos for sad in sads for pos in sad["pos"]]))
    for sad_stat in ["REF", "ALT", "SAD"]:
        sad_h5.create_dataset(sad_stat, data=np.concatenate([sad[sad_stat] for sad in sads], axis=0))

    sad_h5.close()

    # Create vcf containing VCFs from all batches
    if options.vcf:
        vcf_files = natsorted(glob.glob(os.path.join(sad_dir, "chr*/*.vcf")))
        assert len(vcf_files) == options.num_files, f"Expected {options.num_files} sad files for {options.num_files} " \
                                                    f"chromosomes, got {len(vcf_files)} files"

        vcfs = pd.concat([pd.read_table(vcf_file, header=None) for vcf_file in vcf_files])
        vcfs.to_csv(os.path.join(options.out_dir, "all_chrs_snps.vcf"), sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
