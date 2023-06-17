#!/usr/bin/env python
import os
import subprocess
from optparse import OptionParser

import h5py
import numpy as np
import pandas as pd
from collections import defaultdict
from kidney_finemapping.basenji.basenji_utils.dna_io import dna_1hot
from tqdm import tqdm


def main():
    usage = "usage: %prog [options] <motif_file> <ism_summed_pos_shifts_dir> <sad_h5>"
    parser = OptionParser(usage)
    parser.add_option("--min_thresh", dest="min_pval_thresh",
                      default=1e-4, type="float",
                      help="Threshold for min p-value for both ref and alt in fimo_diffs [Default: %default]")
    parser.add_option("--ratio_thresh", dest="ratio_pval_thresh",
                      default=2, type="float",
                      help="Threshold in orders of magnitude for ratio of p-values between alleles [Default: %default]")
    parser.add_option("-t", dest="targets_file",
                      default=None, type="str",
                      help="File specifying target indexes and labels in table format")
    parser.add_option("--stats", dest="stats",
                      default="SUM,MAX", type="str",
                      help="Comma-separated list for heuristic stats: SUM,MAX [Default: %default]")
    parser.add_option("-l", dest="mut_len",
                      default=20, type="int",
                      help="Length of center sequence to mutate used in basenji_ism_summed_pos_shifts [Default: %default]")
    parser.add_option("-o", dest="out_dir",
                      default=None)

    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 3
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    motif_file = args[0]
    ism_summed_pos_shifts_dir = args[1]
    sad_h5_file = args[2]

    os.makedirs(options.out_dir, exist_ok=True)
    options.stats = options.stats.lower().split(",")

    # Read ISM info and SAD h5
    ism_info_file = f"{ism_summed_pos_shifts_dir}/info.h5"
    ism_info = h5py.File(ism_info_file, "r")

    sad_h5 = h5py.File(sad_h5_file, "r")
    sad_rsids = [rsid.decode("utf-8") for rsid in sad_h5["snp"]]
    sad_ref_alleles = [allele.decode("utf-8") for allele in sad_h5["ref_allele"]]
    sad_alt_alleles = [allele.decode("utf-8") for allele in sad_h5["alt_allele"]]
    sad_df = pd.DataFrame({"sad_rsid": sad_rsids, "ref_allele": sad_ref_alleles, "alt_allele": sad_alt_alleles})

    rsids = [rsid.decode("utf-8") for rsid in ism_info["rsid"]]
    ref_seqs = [ref_seq.decode("utf-8") for ref_seq in ism_info["ref_ism_seq"]]
    alt_seqs = [alt_seq.decode("utf-8") for alt_seq in ism_info["alt_ism_seq"]]

    snp_seq_df = pd.DataFrame({"rsid": rsids, "ref_seq": ref_seqs, "alt_seq": alt_seqs})
    snp_seq_df = snp_seq_df.merge(sad_df, how="inner", left_on="rsid", right_on="sad_rsid", validate="1:1").drop("sad_rsid", axis=1)
    snp_seq_df = snp_seq_df.set_index("rsid")

    # Assume symmetric
    mut_up = options.mut_len // 2
    mut_down = options.mut_len // 2

    assert (snp_seq_df["ref_seq"].str[mut_down:-mut_up] == snp_seq_df["ref_allele"]).all(), "ref seqs do not match ref alleles assuming mut length 20"
    assert (snp_seq_df["alt_seq"].str[mut_down:-mut_up] == snp_seq_df["alt_allele"]).all(), "alt seqs do not match alt alleles assuming mut length 20"

    # Create FASTA with reference sequences surrounding each variant
    ref_fasta_file = f"{options.out_dir}/ref_fasta.fa"
    ref_fasta_out = open(ref_fasta_file, "w")

    for i, seq in enumerate(ref_seqs):
        print(f">{rsids[i]}", file=ref_fasta_out)
        print(seq, file=ref_fasta_out)

    ref_fasta_out.close()

    # Create FASTA with alternate sequences surrounding each variant
    alt_fasta_file = f"{options.out_dir}/alt_fasta.fa"
    alt_fasta_out = open(alt_fasta_file, "w")

    for i, seq in enumerate(alt_seqs):
        print(f">{rsids[i]}", file=alt_fasta_out)
        print(seq, file=alt_fasta_out)

    alt_fasta_out.close()

    #######################################################
    # Query ref and alt sequences with FIMO
    ref_fimo_dir = f"{options.out_dir}/ref_fimo_out"
    alt_fimo_dir = f"{options.out_dir}/alt_fimo_out"
    subprocess.call(f"fimo --o {ref_fimo_dir} {motif_file} {ref_fasta_file}", shell=True)
    subprocess.call(f"fimo --o {alt_fimo_dir} {motif_file} {alt_fasta_file}", shell=True)

    ref_fimo_all_file = f"{ref_fimo_dir}/fimo_all.txt"
    ref_fimo_all_out = open(ref_fimo_all_file, "w")
    subprocess.call(f"fimo --thresh 1 --text {motif_file} {ref_fasta_file}",
                    shell=True, stdout=ref_fimo_all_out)

    alt_fimo_all_file = f"{alt_fimo_dir}/fimo_all.txt"
    alt_fimo_all_out = open(alt_fimo_all_file, "w")
    subprocess.call(f"fimo --thresh 1 --text {motif_file} {alt_fasta_file}",
                    shell=True, stdout=alt_fimo_all_out)

    #######################################################
    # Compute metrics for differences between ref and alt FIMO queries
    column_names = ["motif_id", "motif_alt_id", "rsid", "start",
                    "stop", "strand", "score", "p-value", "q-value", "matched_sequence"]

    print("Loading in fimo data...")
    ref_fimo_all = pd.read_table(ref_fimo_all_file, sep="\t", names=column_names, comment="#", keep_default_na=False)
    alt_fimo_all = pd.read_table(alt_fimo_all_file, sep="\t", names=column_names, comment="#", keep_default_na=False)

    # Filter out motif queries that don"t overlap with the variant
    print("Filtering out queries that don"t overlap with the variant...")
    ref_fimo_all = ref_fimo_all.merge(snp_seq_df, how="inner", left_on="rsid", right_index=True, validate="m:1")
    ref_fimo_all["variant_start"] = mut_up + 1  # 1-indexed
    ref_fimo_all["variant_end"] = mut_up + ref_fimo_all["ref_allele"].str.len()  # 1-indexed, inclusive
    overlap = np.maximum(0, np.minimum(ref_fimo_all["variant_end"], ref_fimo_all["stop"]) -
                             np.maximum(ref_fimo_all["start"], ref_fimo_all["variant_start"]) + 1
                         )
    ref_fimo_all = ref_fimo_all[overlap > 0]

    alt_fimo_all = alt_fimo_all.merge(snp_seq_df, how="inner", left_on="rsid", right_index=True, validate="m:1")
    alt_fimo_all["variant_start"] = mut_up + 1  # 1-indexed
    alt_fimo_all["variant_end"] = mut_up + alt_fimo_all["alt_allele"].str.len()  # 1-indexed, inclusive
    overlap = np.maximum(0, np.minimum(alt_fimo_all["variant_end"], alt_fimo_all["stop"]) -
                             np.maximum(alt_fimo_all["start"], alt_fimo_all["variant_start"]) + 1
                         )
    alt_fimo_all = alt_fimo_all[overlap > 0]

    # For each variant, for each motif, get best motif match in the sequence surrounding motif
    ref_fimo_best_matches = ref_fimo_all.groupby(["motif_id", "rsid"]).agg({"p-value": min}).reset_index()
    alt_fimo_best_matches = alt_fimo_all.groupby(["motif_id", "rsid"]).agg({"p-value": min}).reset_index()

    # Compute diff metrics
    print("Computing diff metrics...")
    diff_fimo_all = ref_fimo_best_matches.merge(alt_fimo_best_matches, how="inner",
                                                left_on=["motif_id", "rsid"],
                                                right_on=["motif_id", "rsid"],
                                                suffixes=("_ref", "_alt"))

    # diff_fimo_all["score_diff"] = diff_fimo_all["score_alt"] - diff_fimo_all["score_ref"]
    diff_fimo_all["log10_p-value_ratio"] = np.log10(diff_fimo_all["p-value_alt"]) - np.log10(diff_fimo_all["p-value_ref"])
    diff_fimo_all["min_p-value"] = np.minimum(diff_fimo_all["p-value_alt"], diff_fimo_all["p-value_ref"])
    preference_mapping = pd.Series(data=["REF", "ALT"], index=[0, 1])
    indices = np.argmin([diff_fimo_all["p-value_ref"], diff_fimo_all["p-value_alt"]], axis=0)
    seq_preferences = preference_mapping[indices]
    seq_preferences.index = diff_fimo_all.index
    diff_fimo_all["PREFERENCE"] = seq_preferences

    # Sort and filter by p-value thresholds
    diff_fimo = diff_fimo_all[diff_fimo_all["min_p-value"] < options.min_pval_thresh]
    del diff_fimo_all
    diff_fimo = diff_fimo[np.abs(diff_fimo["log10_p-value_ratio"]) > options.ratio_pval_thresh]
    diff_fimo = diff_fimo.sort_values(by="log10_p-value_ratio", ascending=False, key=np.abs)
    diff_fimo = diff_fimo.reset_index(drop=True)

    # For each motif-SNP pair, get ties for the preferred allele
    ## Define motif_rsid as a convenience column
    print("Retrieving ties...")
    ref_fimo_all["motif_rsid"] = ref_fimo_all["motif_id"] + "_" + ref_fimo_all["rsid"]
    alt_fimo_all["motif_rsid"] = alt_fimo_all["motif_id"] + "_" + alt_fimo_all["rsid"]
    diff_fimo["motif_rsid"] = diff_fimo["motif_id"] + "_" + diff_fimo["rsid"]

    ref_motif_rsids = set(diff_fimo[diff_fimo["PREFERENCE"] == "REF"]["motif_rsid"])
    alt_motif_rsids = set(diff_fimo[diff_fimo["PREFERENCE"] == "ALT"]["motif_rsid"])

    ref_fimo_filtered = ref_fimo_all[ref_fimo_all["motif_rsid"].isin(ref_motif_rsids)]
    del ref_fimo_all

    alt_fimo_filtered = alt_fimo_all[alt_fimo_all["motif_rsid"].isin(alt_motif_rsids)]
    del alt_fimo_all

    ref_mins = ref_fimo_filtered[ref_fimo_filtered["p-value"] == ref_fimo_filtered.groupby("motif_rsid")["p-value"].transform(min)]
    alt_mins = alt_fimo_filtered[alt_fimo_filtered["p-value"] == alt_fimo_filtered.groupby("motif_rsid")["p-value"].transform(min)]

    combined_mins = pd.concat([ref_mins, alt_mins], axis=0)[["motif_alt_id", "start", "stop", "strand", "matched_sequence", "ref_seq", "alt_seq", "motif_rsid"]]
    diff_fimo = diff_fimo.merge(combined_mins, how="inner", left_on="motif_rsid", right_on="motif_rsid", validate="1:m").drop("motif_rsid", axis=1)

    # Rename and reorder columns
    diff_fimo.columns = ["motif_id", "rsid", "best_pval_ref", "best_pval_alt", "log10_pval_ratio",
                         "min_pval", "preference", "motif_alt_id", "start", "stop", "strand",
                         "matched_sequence", "ref_seq", "alt_seq"]

    diff_fimo = diff_fimo[["motif_id", "motif_alt_id", "rsid", "matched_sequence", "ref_seq", "alt_seq", "best_pval_ref", "best_pval_alt",
                           "log10_pval_ratio", "min_pval", "preference", "start", "stop", "strand",
                         ]]

    diff_fimo.to_csv(f"{options.out_dir}/fimo_diffs.txt", index=None, sep="\t", mode="w")

    #######################################################
    # Prioritizing TF-binding motifs based on ISM scores
    print("Prioritizing based on ISM scores...")
    targets_df = pd.read_csv(options.targets_file, sep="\t", index_col=0)
    target_indices = [2, 4, 8, 9, 7, 3, 6, 1, 5, 0]  # Reorder targets

    alleles = ["ref", "alt"]
    motif_scores_dfs = defaultdict(lambda: [])

    for allele in alleles:
        if allele == "ref":
            preds_dir = "ref_ism_preds"
            ism_seq_col = "ref_seq"
        else:
            preds_dir = "alt_ism_preds"
            ism_seq_col = "alt_seq"

        for stat in options.stats:
            motif_scores_df = diff_fimo.copy()

            # Initialize new columns for each target
            for ti in target_indices:
                motif_scores_df[targets_df["identifier"][ti]] = 0

            # Iterate over motifs
            for mi in tqdm(range(diff_fimo.shape[0])):
                rsid = diff_fimo["rsid"].iloc[mi]
                start = diff_fimo["start"].iloc[mi] - 1
                stop = diff_fimo["stop"].iloc[mi]

                # For given rsid, compute scores relative to reference/alternate (original sequence before ISM)
                ism_scores = np.load(f"{ism_summed_pos_shifts_dir}/{preds_dir}/{rsid}.npy")
                seq_1hot = dna_1hot(snp_seq_df.loc[rsid, ism_seq_col])
                ref_scores = np.expand_dims(ism_scores[seq_1hot], axis=1)
                scores_rel_ref = ism_scores - ref_scores

                magnitudes = np.abs(scores_rel_ref)
                if stat == "sum":
                    motif_scores = np.sum(magnitudes[start:stop], axis=(0, 1)) / (stop - start)
                elif stat == "max":
                    motif_scores = np.sum(np.max(magnitudes[start:stop], axis=1), axis=0) / (stop - start)
                else:
                    assert False, "Stat {} not recognized".format(stat)
                for ti, target in enumerate(targets_df["identifier"]):
                    motif_scores_df.loc[mi, f"{target}"] = motif_scores[ti]

            motif_scores_df.insert(diff_fimo.shape[1], "max_score_by_ism",
                                   motif_scores_df.iloc[:, diff_fimo.shape[1]:].max(axis=1))
            motif_scores_df = motif_scores_df.sort_values(by="max_score_by_ism", ascending=False)

            motif_scores_df.to_csv(r"{}".format(
                os.path.join(options.out_dir, "motifs_scored_by_ism_{}_{}.tsv".format(stat, allele))),
                index=None, sep="\t", mode="w")

            motif_scores_dfs[stat].append(motif_scores_df)

    # Merge ref and alt, using only the motif score of the preferred sequence (ref or alt)
    for stat in motif_scores_dfs.keys():
        ref_df, alt_df = motif_scores_dfs[stat]
        ref_df = ref_df[ref_df["preference"] == "REF"]
        alt_df = alt_df[alt_df["preference"] == "ALT"]
        combined_df = pd.concat([ref_df, alt_df])
        combined_df = combined_df.sort_values(by="max_score_by_ism", ascending=False)
        combined_df.to_csv(
            r"{}".format(os.path.join(options.out_dir, "motifs_scored_by_ism_{}_combined.tsv".format(stat))),
            index=None, sep="\t", mode="w")


if __name__ == "__main__":
    main()
