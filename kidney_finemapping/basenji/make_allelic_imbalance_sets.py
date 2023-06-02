#!/usr/bin/env python
import os
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


def main():
    """
    Given allelic imbalance data merged across individuals, create positive and negative sets for analysis.
    """
    usage = "usage: %prog [options] <allelic_imbalance_file> <peaks_file>"
    parser = OptionParser(usage)
    parser.add_option("-o", dest="out_dir",
                      default="variant_sets"),
    parser.add_option("--thresh", dest="thresh",
                      default=0.05, type="float",
                      help="FDR threshold to filter variants by Benjamini-Hochberg [Default: %default]")
    parser.add_option("--neg_mult", dest="neg_mult",
                      default=1, type="int",
                      help="Multiplier for negative set size relative to positive set [Default: %default]")
    parser.add_option("-v", dest="valid_split",
                      default=0.2, type="float",
                      help="Proportion of variants to use for validation set [Default: %default]")
    parser.add_option("--n_bins", dest="n_bins",
                      default=20, type="int",
                      help="Number of intervals for matching log read count distribution")
    parser.add_option("--tubule_peaks", dest="tubule_peaks",
                      default=False, action="store_true",
                      help="Take union of pan tubule peaks [Default: %default]")

    (options, args) = parser.parse_args()

    num_expected_args = 2
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    # Parse arguments
    allelic_imbalance_file = args[0]
    peaks_file = args[1]

    os.makedirs(options.out_dir, exist_ok=True)
    np.random.seed(1)
    #######################################################
    # Get variants of interest
    # Read in allelic imbalance and peaks files
    if options.tubule_peaks:
        # Get union of all pan tubule peaks
        peaks = []
        for cell_type in ["CFH", "PT", "LOH", "DT", "CD"]:
            peak_df = pd.read_table(os.path.join(peaks_file, "{}_peaks.narrowPeak".format(cell_type)), sep="\t",
                                    header=None)
            peaks.append(peak_df)
            peaks = pd.concat(peaks, axis=0)
        peaks = pd.concat(peaks, axis=0)
    else:
        peaks = pd.read_table(peaks_file, sep="\t", header=None)

    ai_df = pd.read_csv(allelic_imbalance_file, sep="\t")

    # Keep only needed columns
    ai_df = ai_df[["Chr", "Pos", "Ref", "Alt", "#Ref", "#Alt", "P-value", "#Ref/(#Ref+#Alt)"]]

    # Filter for total reads > 20
    ai_df["#Total"] = ai_df["#Ref"] + ai_df["#Alt"]
    ai_df = ai_df[ai_df["#Total"] > 20]

    # Filter out indels or triallelic sites
    ai_df = ai_df[~(ai_df["Ref"].str.len() > 1)]
    ai_df = ai_df[~(ai_df["Alt"].str.len() > 1)]

    # Filter for variants within a peak
    select = []
    chr_nums = [chr_num for chr_num in np.unique(ai_df["Chr"])]
    chr_masks = {chr_num: (peaks[0] == "chr{}".format(chr_num)) for chr_num in chr_nums}

    for si in range(ai_df.shape[0]):
        chr_mask = chr_masks[ai_df["Chr"].iloc[si]]
        if np.any((ai_df["Pos"].iloc[si] >= peaks[1][chr_mask]) & (ai_df["Pos"].iloc[si] <= peaks[2][chr_mask])):
            select.append(si)

    ai_df = ai_df.iloc[select]
    ai_df = ai_df.reset_index(drop=True)

    #######################################################
    # Create positive set of variants
    # Compute and filter by FDR on remaining variants by Benjamini-Hochberg
    rejects, qvals, _, _ = sm.stats.multipletests(ai_df["P-value"], alpha=options.thresh, method="fdr_bh")
    ai_df["Q-value"] = qvals
    pos_set = ai_df[rejects]

    #######################################################
    # Create negative set of variants
    # Remove positive set variants and filter by imbalance
    neg_set = pd.concat([ai_df, pos_set], axis=0, ignore_index=True)
    neg_set = neg_set.drop_duplicates(keep=False)
    neg_set = neg_set[neg_set["P-value"] > 0.1]

    # Find a negative set with similar read count distribution to positive set
    pos_log_counts = np.log(pos_set["#Total"])
    neg_log_counts = np.log(neg_set["#Total"])
    bins = np.linspace(np.min(pos_log_counts), np.max(pos_log_counts), num=options.n_bins)
    pos_binned = np.digitize(pos_log_counts, bins=bins)
    neg_binned = np.digitize(neg_log_counts, bins=bins)

    pos_bin_is, pos_counts = np.unique(pos_binned, return_counts=True)
    neg_bin_is, neg_counts = np.unique(neg_binned, return_counts=True)

    assert np.all(np.isin(pos_bin_is, neg_bin_is)), ("No corresponding bin in negative set; "
                                                     "may need to decrease the number of bins")
    match_pos = np.isin(neg_bin_is, pos_bin_is)
    neg_counts = neg_counts[match_pos]

    min_neg_mult = np.min((neg_counts / pos_counts).astype(np.int32))

    assert min_neg_mult >= options.neg_mult, "Neg mult is too large; cannot find enough matching negative variants"

    # Sample bins
    negs = []
    for bin_i in pos_bin_is:
        pos_set_i = pos_set[pos_binned == bin_i]
        neg_set_i = neg_set[neg_binned == bin_i]
        neg_set_i = neg_set_i.sample(n=options.neg_mult * pos_set_i.shape[0],
                                     replace=False, axis=0, random_state=1)
        negs.append(neg_set_i)

    neg_set = pd.concat(negs, axis=0)

    #######################################################
    # Write positive and negative variant sets to output VCF dummy file
    # Prepare variant set VCFs
    pos_set = pos_set.set_index(pos_set["Chr"].astype(str) + ":" + pos_set["Pos"].astype(str))
    pos_set.insert(2, "index", pos_set.index)

    neg_set = neg_set.set_index(neg_set["Chr"].astype(str) + ":" + neg_set["Pos"].astype(str))
    neg_set.insert(2, "index", neg_set.index)

    # Write positive and negative sets to vcfs
    vcf_pos_file = os.path.join(options.out_dir, "pos_set.vcf")
    vcf_out = open(vcf_pos_file, "w")
    print("##fileformat=VCFv4.3", file=vcf_out)
    print("#{}".format("\t".join(pos_set.columns)), file=vcf_out)
    vcf_out.close()
    pos_set.to_csv(r"{}".format(vcf_pos_file), header=None, index=None, sep="\t", mode="a")

    vcf_neg_file = os.path.join(options.out_dir, "neg_set.vcf")
    vcf_out = open(vcf_neg_file, "w")
    print("##fileformat=VCFv4.3", file=vcf_out)
    print("#{}".format("\t".join(neg_set.columns)), file=vcf_out)
    vcf_out.close()
    neg_set.to_csv(r"{}".format(vcf_neg_file), header=None, index=None, sep="\t", mode="a")

    #################################################################
    # Filter VCFs by chromosome
    snps_chrom_dir = os.path.join(options.out_dir, "snps_by_chrom")
    os.makedirs(snps_chrom_dir, exist_ok=True)

    for chrom_num in range(1, 23):
        chrom = "chr" + str(chrom_num)
        pos_vcf_chrom_file = os.path.join(snps_chrom_dir, "pos_set_{}.vcf".format(chrom))
        pos_set[pos_set.iloc[:, 0] == chrom_num].to_csv(pos_vcf_chrom_file, sep="\t", header=False, index=False)

        neg_vcf_chrom_file = os.path.join(snps_chrom_dir, "neg_set_{}.vcf".format(chrom))
        neg_set[neg_set.iloc[:, 0] == chrom_num].to_csv(neg_vcf_chrom_file, sep="\t", header=False, index=False)

    #######################################################
    # Plot read count distributions for positive and negative sets
    plt.title("Read count distributions for variants in positive and negative sets")
    plt.hist(pos_log_counts, alpha=0.5, bins=bins, density=True, label="Positive")
    plt.hist(np.log(neg_set["#Total"]), alpha=0.5, bins=bins, density=True, label="Negative")
    plt.legend()
    plt.savefig(os.path.join(options.out_dir, "read_count_hist.pdf"), dpi=600)
    plt.close("all")

    # All neg vs pos
    plt.title("Read count distributions for variants in positive and all negatives")
    plt.hist(pos_log_counts, alpha=0.5, bins=bins, density=True, label="Positive")
    plt.hist(neg_log_counts, alpha=0.5, bins=bins, density=True, label="All Negative")
    plt.legend()
    plt.savefig(os.path.join(options.out_dir, "read_count_hist_all_neg.pdf"), dpi=600)
    plt.close("all")


if __name__ == "__main__":
    main()
