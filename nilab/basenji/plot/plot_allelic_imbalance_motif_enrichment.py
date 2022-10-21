#!/usr/bin/env python

import os
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    usage = 'usage: %prog [options] <tf_pseudobulk_Pseudobulk_Wilson_TF_analysis csv> <hypergeom_per_motif tsv>'
    parser = OptionParser(usage)
    parser.add_option("--cell_type", dest="cell_type", default=None, help="Cell type for selecting top 300 expressed TFs")
    parser.add_option("-n", dest="n", default=300, help="Number of top expressed TFs to select")
    parser.add_option('-o', dest='out_dir', default='sad_pos_shifts_plots', help='Output directory for script')
    options, args = parser.parse_args()

    # Setup
    num_expected_args = 2
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    tf_analysis_csv = args[0]
    hypergeom_per_motif_tsv = args[1]
    os.makedirs(options.out_dir, exist_ok=True)

    # Read in files
    tf_analysis_df = pd.read_csv(tf_analysis_csv, index_col=0)
    hypergeom_per_motif_df = pd.read_csv(hypergeom_per_motif_tsv, sep="\t", index_col=0)

    hypergeom_per_motif_df = hypergeom_per_motif_df.sort_values(by="hypergeom_qval", ascending=True)  # sort by qval

    # Extract top N tf motifs
    if options.cell_type == "combined":
        tubule_targets = ['CFH', 'PT', 'LOH', 'DT', 'CD']
        for target in tubule_targets:
            raise NotImplementedError("Implement this later")

    top_tfs = set(tf_analysis_df.sort_values(by=options.cell_type, ascending=False)["gene_name"][:options.n].str.lower().values)

    # subset to motifs for top N expressed TFs, taking most significant motif for duplicates (e.g. dealing with var.2)
    motifs = hypergeom_per_motif_df["motif_alt_id"].str.lower().str.replace(r'\([^)]*\)', '', regex=True)  # get rid of (var.2) annotation
    top_n_motifs = set(motifs.values).intersection(top_tfs)
    filtered_hypergeom_df = hypergeom_per_motif_df[hypergeom_per_motif_df["motif_alt_id"].str.lower().isin(top_n_motifs)]
    filtered_hypergeom_df = filtered_hypergeom_df[~filtered_hypergeom_df.duplicated(keep="first")]

    # sanity check: check how many motifs match expressed TFs (not just the top ones)
    motifs_found = set(motifs.values).intersection(tf_analysis_df["gene_name"].str.lower().values)
    print(f"Found {filtered_hypergeom_df.shape[0]} motifs in the top {options.n} motifs")
    print(f"{len(motifs_found)} out of {hypergeom_per_motif_df.shape[0]} motifs were found in TF analysis csv")

    # filter for significant matches by hypergeom for plotting
    filtered_hypergeom_df = filtered_hypergeom_df[filtered_hypergeom_df["hypergeom_qval"] <= 0.05]

    # bar plots
    filtered_hypergeom_df["-log10(hypergeom_qval)"] = -np.log10(filtered_hypergeom_df["hypergeom_qval"])

    # plot all enrichment qvalue
    plt.figure(figsize=(8, 4.8))
    sns.barplot(x="motif_alt_id", y="-log10(hypergeom_qval)", data=filtered_hypergeom_df, color="dimgray")
    # plt.xticks([])
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{options.out_dir}/{options.cell_type}_hypergeom_top_{options.n}.png", dpi=300, bbox_inches="tight")
    plt.show()

    # plot percent positive increased
    plt.figure(figsize=(8, 4.8))
    barplot = sns.barplot(x="motif_alt_id", y="pct_pos_incr", data=filtered_hypergeom_df, color="dimgray")
    plt.axhline(y=0.5, c="orange", ls="--", dashes=(4, 5))
    plt.xticks(rotation=45)

    # add significance asterisks
    filtered_hypergeom_df["binom_sig"] = filtered_hypergeom_df["binom_qval"] <= 0.05
    for p, sig in zip(barplot.patches, filtered_hypergeom_df["binom_sig"].values):
        if sig:
            barplot.text(p.get_x() + p.get_width() / 2., p.get_height(), "*", ha="center")

    plt.tight_layout()
    plt.savefig(f"{options.out_dir}/{options.cell_type}_binom_top_{options.n}.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == '__main__':
    main()
