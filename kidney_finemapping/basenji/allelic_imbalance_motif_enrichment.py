#!/usr/bin/env python
import os
import subprocess
from optparse import OptionParser

import numpy as np
import pandas as pd
import pysam
import statsmodels.api as sm
from scipy.stats import hypergeom, binom_test


def main():
    usage = "usage: %prog [options] <motif_file>"
    parser = OptionParser(usage)
    parser.add_option("--vcf_dir", dest="vcf_dir", default=None, type="str", help="VCF file directory")
    parser.add_option("--ai_dir", dest="ai_dir", default=None, type="str", help="Allelic imbalance file directory")
    parser.add_option("-f", dest="genome_fasta", default="/home/rshuai/research/ni-lab/analysis/genomes/hg38.ml.fa")
    parser.add_option("--min_thresh", dest="min_pval_thresh", default=1e-4, type="float", help="Threshold for min p-value for both ref and alt in fimo_diffs [Default: %default]")
    parser.add_option("--diff_thresh", dest="diff_pval_thresh", default=1, type="float", help="Threshold for magnitude difference in log p-values between ref and alt [Default: %default]")
    parser.add_option("-w", dest="window_size", default=20, type="int", help="Window size centered at SNP for FIMO query")
    parser.add_option("--no_compute_fimo", dest="no_compute_fimo", default=False, action="store_true", help="Compute FIMO before enrichment tests (necessary for first time) [Default: %default]")
    parser.add_option("--no_compute_fimo_diffs", dest="no_compute_fimo_diffs", default=False, action="store_true", help="Compute FIMO diffs before enrichment tests (necessary for first time) [Default: %default]")
    parser.add_option("-o", dest="out_dir", default="ai_motif_enrichment")
    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    motif_file = args[0]

    os.makedirs(options.out_dir, exist_ok=True)
    vcfs = ["pos_set", "neg_set"]

    # Create FASTA with reference sequences surrounding each variant
    num_seqs_per_set = []
    genome_open = pysam.Fastafile(options.genome_fasta)

    left_len = options.window_size // 2
    right_len = options.window_size - left_len

    for vcf in vcfs:
        ref_fimo_dir = os.path.join(options.out_dir, "{}_ref_fimo_out".format(vcf))
        os.makedirs(ref_fimo_dir, exist_ok=True)
        alt_fimo_dir = os.path.join(options.out_dir, "{}_alt_fimo_out".format(vcf))
        os.makedirs(alt_fimo_dir, exist_ok=True)
        ref_fimo_all_file = os.path.join(ref_fimo_dir, "fimo_all.txt")
        alt_fimo_all_file = os.path.join(alt_fimo_dir, "fimo_all.txt")

        if not options.no_compute_fimo:
            ref_fasta_file = os.path.join(options.out_dir, "{}_ref_fasta.fa".format(vcf))
            ref_fasta_out = open(ref_fasta_file, "w")

            alt_fasta_file = os.path.join(options.out_dir, "{}_alt_fasta.fa".format(vcf))
            alt_fasta_out = open(alt_fasta_file, "w")

            vcf_in = open(os.path.join(options.vcf_dir, "{}.vcf".format(vcf)), "r")

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
            subprocess.call("fimo --o {} {} {}"
                            .format(ref_fimo_dir, motif_file, ref_fasta_file),
                            shell=True)

            subprocess.call("fimo --o {} {} {}"
                            .format(alt_fimo_dir, motif_file, alt_fasta_file),
                            shell=True)

            ref_fimo_all_out = open(ref_fimo_all_file, "w")
            subprocess.call("fimo --thresh 1 --text {} {}"
                            .format(motif_file, ref_fasta_file),
                            shell=True, stdout=ref_fimo_all_out)

            alt_fimo_all_out = open(alt_fimo_all_file, "w")
            subprocess.call("fimo --thresh 1 --text {} {}"
                            .format(motif_file, alt_fasta_file),
                            shell=True, stdout=alt_fimo_all_out)

        # Compute metrics for differences between ref and alt FIMO queries
        if not options.no_compute_fimo_diffs:
            column_names = ["motif_id", "motif_alt_id", "sequence_name", "start",
                            "stop", "strand", "score", "p-value", "q-value", "matched_sequence"]

            diff_fimo_all_file = os.path.join(options.out_dir, "{}_fimo_diffs_all.txt".format(vcf))
            diff_fimo_file = os.path.join(options.out_dir, "{}_fimo_diffs.txt".format(vcf))

            # Clear contents
            diff_fimo_all = open(diff_fimo_all_file, "w")
            diff_fimo_all.close()
            diff_fimo = open(diff_fimo_file, "w")
            diff_fimo.close()

            # Read FIMO all files in chunks
            ref_fimo_all = pd.read_table(ref_fimo_all_file, sep="\t", names=column_names, comment="#",
                                         keep_default_na=False, chunksize=1e7)
            alt_fimo_all = pd.read_table(alt_fimo_all_file, sep="\t", names=column_names, comment="#",
                                         keep_default_na=False, chunksize=1e7)
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
    diff_fimo_files = [os.path.join(options.out_dir, "{}_fimo_diffs.txt".format(vcf)) for vcf in vcfs]
    diff_fimos = [pd.read_table(dff, sep="\t", index_col=None, header=0) for dff in diff_fimo_files]
    num_seqs_per_set = [sum(1 for line in open(os.path.join(options.out_dir, "pos_set_ref_fasta.fa"))) // 2,
                        sum(1 for line in open(os.path.join(options.out_dir, "neg_set_ref_fasta.fa"))) // 2]

    # Load allelic imbalance files from file
    ai_files = [os.path.join(options.ai_dir, "{}.vcf".format(vcf)) for vcf in vcfs]
    columns = ["Chr", "Pos", "Seq_name", "Ref", "Alt", "#Ref", "#Alt", "P-value", "#Ref/(#Ref+#Alt)", "#Total",
               "Q-value"]
    ai_dfs = [pd.read_table(ai_file, sep="\t", index_col=None, names=columns, comment="#") for ai_file in ai_files]
    ai_df = pd.concat([ai_dfs[0], ai_dfs[1]], axis=0).reset_index(drop=True)

    # Compute hypergeometric test to find enrichment of disrupted motifs in positive set
    print("========= Performing hypergeometric test ========= ")
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

    pval_df.to_csv(os.path.join(options.out_dir, "hypergeom_per_motif.tsv"), sep="\t")


if __name__ == "__main__":
    main()
