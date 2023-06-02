#!/usr/bin/env python
from __future__ import print_function

import os
from optparse import OptionParser

import pandas as pd


def main():
    """
    Preprocess all finemapped variants for computing SAD scores. Assumes all variants given are in hg38.
    """
    usage = "usage: %prog [options] <finemapped_variants_file>"
    parser = OptionParser(usage)
    parser.add_option("--variant_set", choices=["220513", "220620"], help="Variant set being preprocessed (['220513', '220620'])")
    parser.add_option("-o", dest="out_dir",
                      default="sad",
                      help="Output directory for tables and plots [Default: %default]")
    (options, args) = parser.parse_args()

    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(
            "Incorrect number of arguments, expected {} arguments but got {}".format(num_expected_args, len(args)))
    variants_file = args[0]

    if not os.path.isdir(options.out_dir):
        os.makedirs(options.out_dir)

    #################################################################
    # Filter variants and convert to VCF for SAD predictions
    vcf_file = os.path.join(options.out_dir, "snps.vcf")

    variants_df = pd.read_table(variants_file, sep="\t", header=0)

    if options.variant_set == "220620":
        variants_vcf = variants_df[["hg_38_chrom", "hg38_bp", "SNP", "A1", "A2"]].assign(**{"QUAL": ".", "FILTER": ".", "INFO": "."})
    elif options.variant_set == "220513":
        variants_vcf = variants_df[["chr_hg38", "BP_hg38", "SNP", "A1", "A2"]].assign(**{"QUAL": ".", "FILTER": ".", "INFO": "."})

    # Write header
    vcf_file_out = open(vcf_file, "w")
    print("##fileformat=VCFv4.3", file=vcf_file_out)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_file_out)
    vcf_file_out.close()

    # Write variants
    variants_vcf.to_csv(vcf_file, sep="\t", header=False, index=False, mode="a")

    if options.variant_set == "220620":
        # make hg19 vcf for 220620 variants
        vcf_file_hg19 = os.path.join(options.out_dir, "snps_hg19.vcf")
        variants_df = pd.read_table(variants_file, sep="\t", header=0)

        variants_vcf = variants_df[["hg19_chr", "hg19_bp", "SNP", "A1", "A2"]].assign(**{"QUAL": ".", "FILTER": ".", "INFO": "."})
        variants_vcf["hg19_chr"] = "chr" + variants_vcf["hg19_chr"].astype(str)

        # Write header
        vcf_file_out = open(vcf_file_hg19, "w")
        print("##fileformat=VCFv4.3", file=vcf_file_out)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_file_out)
        vcf_file_out.close()

        # Write variants
        variants_vcf.to_csv(vcf_file_hg19, sep="\t", header=False, index=False, mode="a")

    #################################################################
    # Filter VCF by chromosome
    vcf_df = pd.read_table(vcf_file, sep="\t", comment="#", header=None).iloc[:, 0:5]

    snps_chrom_dir = os.path.join(options.out_dir, "snps_by_chrom")
    os.makedirs(snps_chrom_dir, exist_ok=True)

    for chrom_num in range(1, 23):
        chrom = "chr" + str(chrom_num)
        vcf_chrom_file = os.path.join(snps_chrom_dir, "{}_snps.vcf".format(chrom))
        vcf_df[vcf_df.iloc[:, 0] == chrom].to_csv(vcf_chrom_file, sep="\t", header=False, index=False)

    # X or Y chromosomes
    chroms_XY = {"chrX", "chrY"}
    vcf_chrom_file = os.path.join(snps_chrom_dir, "chrXY_snps.vcf")
    vcf_xy = vcf_df[vcf_df.iloc[:, 0].isin(chroms_XY)]
    if not vcf_xy.empty:
        vcf_xy.to_csv(vcf_chrom_file, sep="\t", header=False, index=False)


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
