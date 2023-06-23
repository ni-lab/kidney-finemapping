#!/usr/bin/env python
import os
import subprocess
from optparse import OptionParser

import numpy as np
import pandas as pd
import pysam
import statsmodels.api as sm
from scipy.stats import binom_test, hypergeom


def main():
    """
    Given positive (allelic imbalance) and negative (non-allelic imbalance) sets of variants and a motif database file,
    compute motif enrichment for each motif in the database.
    """
    usage = "usage: %prog [options] <variant_sets_dir>"
    parser = OptionParser(usage)
    parser.add_option("--motif_file", dest="motif_file",
                      default="./resources/motif_dbs/JASPAR2020_CORE_nonredundant_vertebrates.meme")
    parser.add_option("-f", dest="genome_fasta", default="./resources/genomes/hg38.ml.fa")
    parser.add_option("--min_thresh", dest="min_pval_thresh", default=1e-3, type="float", help="Threshold for min p-value for both ref and alt in fimo_diffs [Default: %default]")
    parser.add_option("--diff_thresh", dest="diff_pval_thresh", default=1, type="float", help="Threshold for magnitude difference in log p-values between ref and alt [Default: %default]")
    parser.add_option("-w", dest="window_size", default=20, type="int", help="Window size centered at SNP for FIMO query")
    parser.add_option("--no_compute_fimo", dest="no_compute_fimo", default=False, action="store_true", help="Do not compute FIMO before enrichment tests (necessary for first time) [Default: %default]")
    parser.add_option("--no_compute_fimo_diffs", dest="no_compute_fimo_diffs", default=False, action="store_true", help="Do not compute FIMO diffs before enrichment tests (must compute for first time) [Default: %default]")
    parser.add_option("-o", dest="out_dir", default="ai_motif_enrichment")
    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    variant_sets_dir = args[0]

    os.makedirs(options.out_dir, exist_ok=True)

    # Create FASTA with reference sequences surrounding each variant
    vcfs = ["pos_set", "neg_set"]
    num_seqs_per_set = []
    genome_open = pysam.Fastafile(options.genome_fasta)

    left_len = options.window_size // 2
    right_len = options.window_size - left_len

    for vcf in vcfs:
        ref_fimo_dir = f"{options.out_dir}/{vcf}_ref_fimo_out"
        os.makedirs(ref_fimo_dir, exist_ok=True)
        ref_fimo_all_file = f"{ref_fimo_dir}/fimo_all.txt"

        alt_fimo_dir = f"{options.out_dir}/{vcf}_alt_fimo_out"
        os.makedirs(alt_fimo_dir, exist_ok=True)
        alt_fimo_all_file = f"{alt_fimo_dir}/fimo_all.txt"

        if not options.no_compute_fimo:
            ref_fasta_file = f"{options.out_dir}/{vcf}_ref_fasta.fa"
            ref_fasta_out = open(ref_fasta_file, "w")

            alt_fasta_file = f"{options.out_dir}/{vcf}_alt_fasta.fa"
            alt_fasta_out = open(alt_fasta_file, "w")

            vcf_in = open(f"{variant_sets_dir}/{vcf}.vcf", "r")

            line = vcf_in.readline()
            while line[0] == "#":
                line = vcf_in.readline()

            num_seqs = 0
            while line:
                chr_num, pos, label, ref_allele, alt_allele, *_ = line.split("\t")
                assert len(ref_allele) == 1, "Ref allele must have length 1"
                assert len(alt_allele) == 1, "Alt allele must have length 1"

                chrom = "chr" + chr_num
                seq_ref = genome_open.fetch(chrom, int(pos) - 1 - left_len, int(pos) - 1 + right_len)

                assert seq_ref[left_len] == ref_allele, \
                    "Reference allele does not match reference genome for SNP {}".format(label)

                print(">{}".format(label), file=ref_fasta_out)
                print(seq_ref, file=ref_fasta_out)

                seq_alt = seq_ref[:left_len] + alt_allele + seq_ref[left_len + 1:]
                print(">{}".format(label), file=alt_fasta_out)
                print(seq_alt, file=alt_fasta_out)

                num_seqs += 1
                line = vcf_in.readline()

            num_seqs_per_set.append(num_seqs)

            ref_fasta_out.close()
            alt_fasta_out.close()

            vcf_in.close()

            # Query ref and alt sequences with FIMO
            subprocess.call(f"fimo --o {ref_fimo_dir} {options.motif_file} {ref_fasta_file}", shell=True)
            subprocess.call(f"fimo --o {alt_fimo_dir} {options.motif_file} {alt_fasta_file}", shell=True)

            ref_fimo_all_out = open(ref_fimo_all_file, "w")
            subprocess.call(f"fimo --thresh 1 --text {options.motif_file} {ref_fasta_file}", shell=True, stdout=ref_fimo_all_out)
            ref_fimo_all_out.close()

            alt_fimo_all_out = open(alt_fimo_all_file, "w")
            subprocess.call(f"fimo --thresh 1 --text {options.motif_file} {alt_fasta_file}", shell=True, stdout=alt_fimo_all_out)
            alt_fimo_all_out.close()

        # Compute metrics for differences between ref and alt FIMO queries
        if not options.no_compute_fimo_diffs:
            diff_fimo_all_file = f"{options.out_dir}/{vcf}_fimo_diffs_all.txt"
            diff_fimo_file = f"{options.out_dir}/{vcf}_fimo_diffs.txt"

            # Clear contents
            diff_fimo_all = open(diff_fimo_all_file, "w")
            diff_fimo_all.close()
            diff_fimo = open(diff_fimo_file, "w")
            diff_fimo.close()

            # Read FIMO all files in chunks
            ref_fimo_all = pd.read_csv(ref_fimo_all_file, sep="\t", header=0, keep_default_na=False, chunksize=1e7)
            alt_fimo_all = pd.read_csv(alt_fimo_all_file, sep="\t", header=0, keep_default_na=False, chunksize=1e7)
            append_header = True
            for ref_chunk, alt_chunk in zip(ref_fimo_all, alt_fimo_all):
                assert len(ref_chunk) == len(
                    alt_chunk), "Different number of queries in ref and alt fimo for some reason..."

                diff_chunk = ref_chunk.copy()
                diff_chunk["score_diff"] = alt_chunk["score"] - ref_chunk["score"]
                diff_chunk["log10_p-value_ratio"] = np.log10(ref_chunk["p-value"]) - np.log10(alt_chunk["p-value"])
                diff_chunk["min_p-value"] = np.minimum(ref_chunk["p-value"], alt_chunk["p-value"])
                preference_mapping = pd.Series(data=["REF", "ALT"], index=[0, 1])
                indices = np.argmin([ref_chunk["p-value"], alt_chunk["p-value"]], axis=0)
                seq_preferences = preference_mapping[indices]
                seq_preferences.index = diff_chunk.index
                diff_chunk["PREFERENCE"] = seq_preferences

                # Write all to file
                diff_chunk = diff_chunk.drop(["score", "p-value", "q-value"], axis=1)
                diff_chunk.to_csv(diff_fimo_all_file, index=None, sep="\t", mode="a", header=append_header)

                # Filter by p-value threshold
                diff_chunk = diff_chunk[diff_chunk["min_p-value"] < options.min_pval_thresh].reset_index(drop=True)

                # Filter by magnitude of log ratio
                diff_chunk["abs_pdiff"] = np.abs(diff_chunk["log10_p-value_ratio"])
                diff_chunk = diff_chunk[diff_chunk["abs_pdiff"] > options.diff_pval_thresh]

                # Save to csv
                diff_chunk.to_csv(diff_fimo_file, index=None, sep="\t", mode="a", header=append_header)

                append_header = False

    genome_open.close()

    # Load filtered diff fimos from file
    diff_fimo_files = [f"{options.out_dir}/{vcf}_fimo_diffs.txt" for vcf in vcfs]
    diff_fimos = [pd.read_table(dff, sep="\t", index_col=None, header=0) for dff in diff_fimo_files]
    num_seqs_per_set = [sum(1 for _ in open(f"{options.out_dir}/pos_set_ref_fasta.fa")) // 2,
                        sum(1 for _ in open(f"{options.out_dir}/neg_set_ref_fasta.fa")) // 2]


    # Load allelic imbalance files from file
    ai_files = [f"{variant_sets_dir}/{vcf}.vcf" for vcf in vcfs]
    columns = ["Chr", "Pos", "Seq_name", "Ref", "Alt", "#Ref", "#Alt", "P-value", "#Ref/(#Ref+#Alt)", "#Total",
               "Q-value"]
    ai_dfs = [pd.read_table(ai_file, sep="\t", index_col=None, names=columns, comment="#") for ai_file in ai_files]
    ai_df = pd.concat([ai_dfs[0], ai_dfs[1]], axis=0).reset_index(drop=True)

    # Compute hypergeometric test to find enrichment of disrupted motifs in positive set
    print("========= Performing hypergeometric tests ========= ")
    pos_diff_fimo, neg_diff_fimo = diff_fimos

    num_pos_seqs, num_neg_seqs = num_seqs_per_set
    total_seqs = num_pos_seqs + num_neg_seqs
    n_pos_matches, n_neg_matches = pos_diff_fimo.shape[0], neg_diff_fimo.shape[0]
    total_matches = n_pos_matches + n_neg_matches

    # Compute hypergeometric test separately for each motif
    pos_diff_fimo["set"] = "POS"
    neg_diff_fimo["set"] = "NEG"
    test_fimo = pd.concat([pos_diff_fimo, neg_diff_fimo], axis=0).reset_index(drop=True)

    # Merge FIMO diffs with allelic imbalance data
    test_fimo = test_fimo.merge(ai_df, left_on="sequence_name", right_on="Seq_name")

    # Take the most disrupted match per sequence per motif to ensure only one match per sequence
    test_fimo = test_fimo.iloc[test_fimo.groupby(["sequence_name", "motif_id"])["abs_pdiff"].idxmax()]

    motif_to_pval = {}
    for (motif, motif_alt), matches in test_fimo.groupby(["motif_id", "motif_alt_id"]):
        # Hypergeometric test
        pos_matches = matches[matches["set"] == "POS"]
        neg_matches = matches[matches["set"] == "NEG"]
        k, M, n, N = pos_matches.shape[0], total_seqs, num_pos_seqs, matches.shape[0]
        hypergeom_pval = hypergeom.sf(k - 1, M, n, N)
        pct_pos_matched = pos_matches.shape[0] / num_pos_seqs
        pct_neg_matched = neg_matches.shape[0] / num_neg_seqs

        # Binomial test within positive set (if we have matches in the positive set)
        # 1 if sequence with motif match has higher read count, 0 if lower.
        if pos_matches.shape[0] > 0:
            pos_motif_dirs = np.where(pos_matches["PREFERENCE"] == "REF", pos_matches["#Ref/(#Ref+#Alt)"] > 0.5,
                                      pos_matches["#Ref/(#Ref+#Alt)"] < 0.5)
            num_pos_incr = np.sum(pos_motif_dirs)
            binom_pval = binom_test(num_pos_incr, pos_matches.shape[0], p=0.5, alternative="two-sided")
            pct_pos_incr = num_pos_incr / pos_matches.shape[0]
            pct_pos_decr = 1 - pct_pos_incr
        else:
            pct_pos_incr, pct_pos_decr, binom_pval = None, None, 1
        motif_to_pval[motif] = [motif_alt, num_pos_seqs, pct_pos_matched, num_neg_seqs, pct_neg_matched,
                                pos_matches.shape[0], pct_pos_incr, pct_pos_decr, hypergeom_pval, binom_pval]

    pval_df = pd.DataFrame.from_dict(motif_to_pval, orient="index",
                                     columns=["motif_alt_id", "num_pos_seqs", "pct_pos_matched", "num_neg_seqs",
                                              "pct_neg_matched",
                                              "num_pos_matches", "pct_pos_incr", "pct_pos_decr", "hypergeom_pval",
                                              "binom_pval"])
    pval_df["hypergeom_qval"] = sm.stats.multipletests(pval_df["hypergeom_pval"], alpha=0.05, method="fdr_bh")[1]
    pval_df["binom_qval"] = sm.stats.multipletests(pval_df["binom_pval"], alpha=0.05, method="fdr_bh")[1]
    pval_df = pval_df.sort_values("hypergeom_pval")

    pval_df.to_csv(f"{options.out_dir}/hypergeom_per_motif.tsv", sep="\t")


if __name__ == "__main__":
    main()
