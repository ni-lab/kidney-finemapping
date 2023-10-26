#!/usr/bin/env python

from __future__ import print_function

import glob
import os
from optparse import OptionParser
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from statannot import add_stat_annotation
from tqdm import tqdm


################################################################################
# main
################################################################################
def main():
    """
    Given a list of finemapped variants, their SAD scores, and allelic imbalance data, plot various statistics binned by PIP.
    """
    usage = 'usage: %prog [options] <finemapped_variants_tsv> <sad_h5> <imbalance_data_dir>'
    parser = OptionParser(usage)
    parser.add_option("--pip_bins", dest="pip_bins", default="0,0.05,0.3,0.7,1.0", help="comma-separated list of PIP bins")
    parser.add_option('-t', dest='targets_file',
                      default=None,
                      type=str,
                      help='File specifying target indexes and labels in table format')
    parser.add_option('-o', dest='out_dir',
                      default=None,
                      help='Output directory for plots')
    (options, args) = parser.parse_args()

    num_expected_args = 3
    if len(args) != num_expected_args:
        parser.error(
            "Incorrect number of arguments, expected {} arguments but got {}".format(num_expected_args, len(args)))

    # Parse arguments
    finemapped_variants_tsv = args[0]
    sad_h5_file = args[1]
    allelic_imbalance_dir = args[2]
    Path(options.out_dir).mkdir(parents=True, exist_ok=True)
    np.random.seed(0)


    # Read in dataframe and parse pip bins
    targets_df = pd.read_csv(options.targets_file, sep='\t', index_col=0)
    pip_bins = list(map(float, options.pip_bins.split(",")))

    # Read PIP data and assign PIP bins
    df_pip = pd.read_csv(finemapped_variants_tsv, sep="\t")
    df_pip["pip_category"] = df_pip["susie10_noprior_PIP"].apply(lambda x: assign_pip_category(x, pip_bins))

    assert ~np.any(df_pip.duplicated(subset=["SNP", "A1", "A2"])), "Expecting no duplicate SNPs in PIP df"
    pip_categories = sorted(np.unique(df_pip["pip_category"]))

    ### Plot allelic imbalance data ###
    imbalance_targets = [os.path.basename(x).replace("all_", "").replace("q10.tsv", "") for x in
                         glob.glob(f"{allelic_imbalance_dir}/*.tsv")]
    ai_out_dir = f"{options.out_dir}/allelic_imbalance"

    Path(ai_out_dir).mkdir(parents=True, exist_ok=True)
    df_imb = df_pip.copy()

    # Add SAD scores to df, making sure to change SAD scores to proper direction
    # for convenience, set index to snp with both alleles sorted for easier indexing
    df_pip.index = df_pip["SNP"] + "_" + (df_pip["A1"] + "_" + df_pip["A2"]).str.split("_").apply(sorted).str.join("_")
    sad_h5 = h5py.File(sad_h5_file, "r")

    snps = [snp.decode("utf-8") for snp in np.array(sad_h5["snp"])]
    refs = [ref.decode("utf-8") for ref in np.array(sad_h5["ref_allele"])]
    alts = [alt.decode("utf-8") for alt in np.array(sad_h5["alt_allele"])]
    sads = [sad for sad in np.array(sad_h5["SAD"])]
    ref_preds = [sad for sad in np.array(sad_h5["REF"])]
    alt_preds = [sad for sad in np.array(sad_h5["ALT"])]

    data = {"idx": [], "sad": [], "ref_pred": [], "alt_pred": []}
    for si in tqdm(range(sad_h5["snp"].shape[0])):
        snp = snps[si]
        ref, alt = refs[si], alts[si]
        sad, ref_pred, alt_pred = sads[si], ref_preds[si].squeeze(), alt_preds[si].squeeze()

        idx = snp + "_" + "_".join(sorted([ref, alt]))
        df_pip_row = df_pip.loc[idx]

        if isinstance(df_pip_row, pd.DataFrame):
            # if we find two matches, use the one with the correct alleles and change index of original df
            assert df_pip_row.shape[0] == 2, "Expecting only two duplicate indices"

            found_correct_alleles = False
            for i, (_, row) in enumerate(df_pip_row.iterrows()):
                if row.loc["A1"] == ref and row.loc["A2"] == alt:
                    data["sad"].append(sad)
                    data["ref_pred"].append(ref_pred)
                    data["alt_pred"].append(alt_pred)

                    idx_new = snp + "_" + ref + "_" + alt
                    data["idx"].append(idx_new)

                    # rename index based on which row we used
                    df_pip_index = df_pip.index.tolist()
                    df_pip_index[np.where(df_pip.index == idx)[0][i]] = idx_new
                    df_pip.index = df_pip_index
                    found_correct_alleles = True
            assert found_correct_alleles, f"Did not find correct alleles for duplicate snp {snp}"
            continue

        data["idx"].append(idx)
        if df_pip_row.loc["A1"] == ref and df_pip_row.loc["A2"] == alt:
            data["sad"].append(sad)
            data["ref_pred"].append(ref_pred)
            data["alt_pred"].append(alt_pred)
        elif df_pip_row.loc["A2"] == ref and df_pip_row.loc["A1"] == alt:
            data["sad"].append(sad)
            data["ref_pred"].append(ref_pred)
            data["alt_pred"].append(alt_pred)
        else:
            assert False, f"SNP {idx} in sad_h5 does not have matching ref and alt alleles to finemapped variants df"

    sad_df = pd.DataFrame(data).set_index("idx")
    df_pip[[f"SAD_{id}" for id in targets_df["identifier"]]] = pd.DataFrame(sad_df["sad"].tolist(), index=sad_df.index)
    df_pip[[f"REF_{id}" for id in targets_df["identifier"]]] = pd.DataFrame(sad_df["ref_pred"].tolist(), index=sad_df.index)
    df_pip[[f"ALT_{id}" for id in targets_df["identifier"]]] = pd.DataFrame(sad_df["alt_pred"].tolist(), index=sad_df.index)

    # PLOTS
    df_pip['pip_category_categorical'] = df_pip['pip_category'].astype(pd.CategoricalDtype(categories=pip_categories, ordered=True))  # deal with statannot bug by making sure PIP categories in df are ordered
    df_pip = df_pip.sort_values(by='pip_category_categorical')

    # # Box plot by SAD score
    # for id in targets_df["identifier"]:
    #     y = df_pip[f"SAD_{id}"]
    #     x = df_pip["pip_category"]
    #     outlier_mask = (y > y.quantile(q=0.02)) & (y < y.quantile(q=0.99))
    #     y_plot = y[outlier_mask]
    #     x_plot = x[outlier_mask]
    #
    #     plt.figure()
    #     ax = sns.boxplot(x=x_plot, y=y_plot, order=pip_categories)
    #     add_stat_annotation(ax, x=x, y=y,
    #                         box_pairs=[(pip_categories[0], x) for x in pip_categories if x != pip_categories[0]],
    #                         test="Mann-Whitney", text_format="full",
    #                         loc="outside")
    #     plt.ylabel(f"SAD_{id}")
    #     plt.tight_layout()
    #     plt.savefig(f"{options.out_dir}/pip_vs_abs_{id}.pdf", dpi=300, bbox_inches="tight")
    #     plt.show()

    # Box plot by |SAD score|
    for id in targets_df["identifier"]:
        y = np.abs(df_pip[f"SAD_{id}"])
        x = df_pip["pip_category"]
        outlier_mask = (y > y.quantile(q=0.02)) & (y < y.quantile(q=0.99))
        y_plot = y[outlier_mask]
        x_plot = x[outlier_mask]

        plt.figure()
        ax = sns.boxplot(x=x_plot, y=y_plot, order=pip_categories)
        add_stat_annotation(ax, x=x, y=y,
                            box_pairs=[(pip_categories[0], x) for x in pip_categories if x != pip_categories[0]],
                            test="Mann-Whitney-ls", text_format="full",
                            loc="outside")
        plt.ylabel(f"|SAD_{id}|")
        plt.tight_layout()
        plt.savefig(f"{options.out_dir}/pip_vs_abs_sad_{id}.pdf", dpi=300, bbox_inches="tight")
        plt.show()

    # Box plot by max |SAD score| across tubule
    tubule_targets = ['CFH', 'PT', 'LOH', 'DT', 'CD']
    y = np.abs(df_pip[[f"SAD_{target}" for target in tubule_targets]]).max(axis=1)
    x = df_pip["pip_category"]
    outlier_mask = (y > y.quantile(q=0.02)) & (y < y.quantile(q=0.99))
    y_plot = y[outlier_mask]
    x_plot = x[outlier_mask]

    plt.figure()
    ax = sns.boxplot(x=x_plot, y=y_plot, order=pip_categories)
    add_stat_annotation(ax, x=x, y=y,
                        box_pairs=[(pip_categories[0], x) for x in pip_categories if x != pip_categories[0]],
                        test="Mann-Whitney-ls", text_format="full",
                        loc="outside")
    plt.ylabel(f"max |SAD| across [CFH, PT, LOH, DT, CD]")
    plt.tight_layout()
    plt.savefig(f"{options.out_dir}/pip_vs_abs_sad_tubule.pdf", dpi=300, bbox_inches="tight")
    plt.show()


    # Box plot by predicted allelic imbalance
    for id in targets_df["identifier"]:
        y = np.abs(0.5 - (df_pip[f"REF_{id}"] / (df_pip[f"REF_{id}"] + df_pip[f"ALT_{id}"])))
        x = df_pip["pip_category"]
        outlier_mask = (y > y.quantile(q=0.02)) & (y < y.quantile(q=0.99))
        y_plot = y[outlier_mask]
        x_plot = x[outlier_mask]

        plt.figure()
        ax = sns.boxplot(x=x_plot, y=y_plot, order=pip_categories)
        add_stat_annotation(ax, x=x, y=y,
                            box_pairs=[(pip_categories[0], x) for x in pip_categories if x != pip_categories[0]],
                            test="Mann-Whitney-ls", text_format="full",
                            loc="outside")
        plt.ylabel(f"predicted |0.5 - IMB_{id}|")
        plt.tight_layout()
        plt.savefig(f"{options.out_dir}/pip_vs_abs_pred_imb_{id}.pdf", dpi=300, bbox_inches="tight")
        plt.show()

    # Box plot by max predicted allellic imbalance across tubule
    tubule_targets = ['CFH', 'PT', 'LOH', 'DT', 'CD']
    ref_scores = df_pip[[f"REF_{target}" for target in tubule_targets]]
    alt_scores = df_pip[[f"ALT_{target}" for target in tubule_targets]]
    y = pd.Series(np.abs(0.5 - (ref_scores.values / (ref_scores.values + alt_scores.values))).max(axis=1), index=ref_scores.index)

    x = df_pip["pip_category"]
    outlier_mask = (y > y.quantile(q=0.02)) & (y < y.quantile(q=0.99))
    y_plot = y[outlier_mask]
    x_plot = x[outlier_mask]

    plt.figure()
    ax = sns.boxplot(x=x_plot, y=y_plot, order=pip_categories)
    add_stat_annotation(ax, x=x, y=y,
                        box_pairs=[(pip_categories[0], x) for x in pip_categories if x != pip_categories[0]],
                        test="Mann-Whitney-ls", text_format="full",
                        loc="outside")
    plt.ylabel(f"max predicted |0.5 - IMB| across [CFH, PT, LOH, DT, CD]")
    plt.tight_layout()
    plt.savefig(f"{options.out_dir}/pip_vs_abs_pred_imb_tubule.pdf", dpi=300, bbox_inches="tight")
    plt.show()

    # Save number of variants of each type
    df_pip.groupby("pip_category").count()[["SAD_PT"]].to_csv(f"{options.out_dir}/variant_counts.tsv", sep="\t", header=True, index=True)


def assign_pip_category(pip, pip_bins):
    for i in range(len(pip_bins) - 2):
        if pip_bins[i] <= pip < pip_bins[i + 1]:
            return f"[{pip_bins[i]}, {pip_bins[i + 1]})"

    # last bin is both endpoints inclusive
    if pip_bins[-2] <= pip <= pip_bins[-1]:
        return f"[{pip_bins[-2]}, {pip_bins[-1]}]"
    assert False, f"pip score of {pip} outside of pip bins {pip_bins}"


if __name__ == '__main__':
    main()
