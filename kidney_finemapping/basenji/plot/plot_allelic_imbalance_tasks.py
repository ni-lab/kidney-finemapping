# !/usr/bin/env python
import bisect
import os
from collections import defaultdict
from optparse import OptionParser
from pathlib import Path
from typing import Optional

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve, roc_curve

from kidney_finemapping.basenji.plot.compare_auc_delong_xu import \
    delong_roc_variance
from kidney_finemapping.kidney_utils import *


def main():
    """
    Given allelic imbalance variant sets and SAD files, plot ROCs and PR curves for direction and classification.
    - args[0] <variant_set_dir> - Directory containing variant sets and SAD files
    """
    usage = "usage: %prog [options] <variant_set_dir>"
    parser = OptionParser(usage)
    parser.add_option("-o", dest="out_dir",
                      default="plots"),
    parser.add_option("-t", dest="targets_file",
                      default=None, type="str",
                      help="File specifying target indexes and labels in table format")

    (options, args) = parser.parse_args()
    num_expected_args = 1

    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    # Parse arguments
    variant_sets_dir = args[0]

    # Setup
    os.makedirs(options.out_dir, exist_ok=True)

    # Parse variant set directory name
    target, neg_mult, thresh = Path(variant_sets_dir).name.split("_")
    neg_mult = neg_mult.replace("neg", "").replace("x", "")
    thresh = thresh.replace("q", "")

    # Load files
    col_names = ["chr", "pos", "index", "ref", "alt", "ref_reads", "alt_reads", "Pvalue", "imbalance", "total_reads", "Qvalue"]
    pos_vcf_df = pd.read_table(f"{variant_sets_dir}/pos_sad/all_chrs/all_chrs_snps.vcf", sep="\t", names=col_names, index_col=None)
    neg_vcf_df = pd.read_table(f"{variant_sets_dir}/neg_sad/all_chrs/all_chrs_snps.vcf", sep="\t", names=col_names, index_col=None)

    pos_sad = h5py.File(f"{variant_sets_dir}/pos_sad/all_chrs/sad.h5", mode="r")
    neg_sad = h5py.File(f"{variant_sets_dir}/neg_sad/all_chrs/sad.h5", mode="r")

    targets_df = pd.read_csv(options.targets_file, sep="\t", index_col=0)

    # Plot read count distributions for positive and negative sets
    plot_read_count(pos_vcf_df, neg_vcf_df, target, neg_mult, thresh, options.out_dir)

    # Plot allelic imbalance direction with cross cell type predictions
    metrics_df = pd.DataFrame(columns=targets_df["identifier"][PLOT_KIDNEY_TARGET_INDICES])  # initialize metrics df to store AUC / AUC std
    metrics_df = plot_roc(pos_vcf_df, pos_sad, targets_df, target, metrics_df, options.out_dir, mode="direction")

    # Plot allelic imbalance classification with cross cell type predictions
    metrics_df = plot_roc(pos_vcf_df, pos_sad, targets_df, target, metrics_df, options.out_dir, mode="class",
                          neg_vcf_df=neg_vcf_df, neg_sad=neg_sad)
    metrics_df = plot_prc_class(pos_vcf_df, pos_sad,
                                neg_vcf_df, neg_sad,
                                targets_df, target, metrics_df, options.out_dir)

    # Save standard error and AUC metrics to file
    metrics_df.to_csv(os.path.join(options.out_dir, "metrics.tsv"), sep="\t", index=True, header=True)


def plot_read_count(
    pos_vcf_df: pd.DataFrame,
    neg_vcf_df: pd.DataFrame,
    target: str,
    neg_mult: str,
    thresh: str,
    out_dir: str
) -> None:
    """
    Plot read count distributions for positive and negative sets.

    Args:
        pos_vcf_df: DataFrame containing positive variant set data.
        neg_vcf_df: DataFrame containing negative variant set data.
        target: The target identifier.
        neg_mult: Negative multiplication factor for scaling.
        thresh: The threshold for filtering.
        out_dir: The directory where to save the plot.

    Returns:
        None
    """
    plt.figure()
    plt.title(f"Read counts for pos and neg sets for {target} (neg_mult={neg_mult}, thresh={thresh})")
    x_max = np.max(pos_vcf_df["total_reads"]) + 1
    plt.hist(pos_vcf_df["total_reads"], alpha=0.5, bins=np.arange(0, x_max, (x_max // 20) + 1), density=True, label="Positive set")
    plt.hist(neg_vcf_df["total_reads"], alpha=0.5, bins=np.arange(0, x_max, (x_max // 20) + 1), density=True, label="Negative set")
    plt.legend()
    plt.savefig(f"{out_dir}/read_counts_{target}_neg{neg_mult}x_q{thresh}.pdf", dpi=300)
    plt.close("all")



def plot_roc(
    pos_vcf_df: pd.DataFrame,
    pos_sad: h5py.File,
    targets_df: pd.DataFrame,
    target: str,
    metrics_df: pd.DataFrame,
    out_dir: str,
    mode: str,
    neg_vcf_df: Optional[pd.DataFrame] = None,
    neg_sad: Optional[h5py.File] = None,
) -> pd.DataFrame:
    """
    Plot ROC curve either for predicting allelic imbalance direction or classifying allelic imbalance from non-allelic imbalance set.
    Updates metrics_df with AUROC and standard deviation values.

    Args:
        pos_vcf_df: DataFrame containing variant set data.
        pos_sad: File with SAD values for each variant.
        targets_df: DataFrame containing target indexes and labels.
        target: The target identifier.
        metrics_df: Dataframe with AUROC and standard deviation values for each target.
        out_dir: The directory where to save the plot.
        mode: one of ["direction", "class"] to indicate whether to plot direction or classification ROC.
        neg_vcf_df: if mode is "class", then this is the non-allelic imbalance set DataFrame.
        neg_sad: if mode is "class", then this is the non-allelic imbalance set SAD file.

    Returns:
        metrics_df: DataFrame with AUROC and standard deviation values for each target.
    """
    plt.figure(figsize=(5.5, 5))

    # Get labels
    if mode == "direction":
        imb_labels = pos_vcf_df["imbalance"] > 0.5
    elif mode == "class":
        imb_labels = np.concatenate([np.ones(pos_vcf_df.shape[0]),
                                     np.zeros(neg_vcf_df.shape[0])])
    else:
        f"Invalid mode: {mode}"

    aurocs = []  # AUROC for direction for sorting legend
    stds = []  # std of auroc for direction
    fpr_thresh = 0.05

    for ti in PLOT_KIDNEY_TARGET_INDICES:
        target_i = targets_df["identifier"][ti]

        if mode == "direction":
            imb_preds = get_imb_from_sad(pos_sad, ti)
        elif mode == "class":
            pos_imbs = get_imb_from_sad(pos_sad, ti)
            neg_imbs = get_imb_from_sad(neg_sad, ti)

            imb_preds = np.concatenate([pos_imbs, neg_imbs])
            imb_preds = np.abs(imb_preds - 0.5)  # classify symmetrically around 0.5

        fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
        roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
        aurocs.append(roc_auc)
        stds.append(auc_cov)

        if target_i in NON_TUBULE_EPITHELIAL_TARGETS:
            color = "lightgray"
        elif target_i in ROC_COLORS:
            color = ROC_COLORS[target_i]
        else:
            assert False, f"Target {target_i} not in ROC_COLORS"

        legend_label = targets_df["identifier"][ti]

        # rename CFH to PEC
        if legend_label == "CFH":
            legend_label = "PEC"

        if target_i == target:
            # plot main line thicker
            plt.plot(fpr, tpr, color=color,
                        lw=3, label="{} AUROC={:.3f}".format(legend_label, roc_auc))

            # get tpr and allelic imbalance threshold when fpr is 0.05
            th_i = bisect.bisect(fpr, fpr_thresh) - 1
            tpr_th = tpr[th_i]
            ai_th = thresholds[th_i]
        else:
            plt.plot(fpr, tpr, color=color,
                        lw=1, label="{} AUROC={:.3f}".format(legend_label, roc_auc))

    plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC, predicting allelic imbalance {mode} \n (FPR: {fpr_thresh}, TPR: {tpr_th:.3f}, AI: {ai_th:.3f})")

    # Reorder by auc
    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.argsort(aurocs)[::-1]
    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])

    # Save plot
    plt.savefig(f"{out_dir}/roc_{mode}_cross.pdf", dpi=300)
    plt.close("all")
    metrics_df.loc[f"{target}_auroc_{mode}"] = aurocs
    metrics_df.loc[f"{target}_std_{mode}"] = stds

    return metrics_df


def plot_prc_class(
    pos_vcf_df: pd.DataFrame,
    pos_sad: h5py.File,
    neg_vcf_df: pd.DataFrame,
    neg_sad: h5py.File,
    targets_df: pd.DataFrame,
    target: str,
    metrics_df: pd.DataFrame,
    out_dir: str,
) -> pd.DataFrame:
    """
    Plot Precision-Recall Curve for predicting allelic imbalance class.
    Updates metrics_df with AUPRC values.

    Args:
        pos_vcf_df: DataFrame containing variant set data.
        pos_sad: File with SAD values for each variant.
        neg_vcf_df: The non-allelic imbalance set DataFrame.
        neg_sad: The non-allelic imbalance set SAD file.
        targets_df: DataFrame containing target indexes and labels.
        target: The target identifier.
        metrics_df: Dataframe with AUC values for each target.
        out_dir: The directory where to save the plot.

    Returns:
        metrics_df: DataFrame with AUPRC values for each target.
    """
    plt.figure(figsize=(5.5, 5))

    # Get labels
    imb_labels = np.concatenate([np.ones(pos_vcf_df.shape[0]),
                                 np.zeros(neg_vcf_df.shape[0])])

    auprcs = []  # AUPRC for direction for sorting legend
    num_pos = np.sum(imb_labels)
    num_total = len(imb_labels)

    for ti in PLOT_KIDNEY_TARGET_INDICES:
        target_i = targets_df["identifier"][ti]

        pos_imbs = get_imb_from_sad(pos_sad, ti)
        neg_imbs = get_imb_from_sad(neg_sad, ti)
        imb_preds = np.concatenate([pos_imbs, neg_imbs])
        imb_preds = np.abs(imb_preds - 0.5)  # classify symmetrically around 0.5

        precision, recall, thresholds = precision_recall_curve(imb_labels, imb_preds)
        auprc = auc(recall, precision)
        auprcs.append(auprc)

        if target_i in NON_TUBULE_EPITHELIAL_TARGETS:
            color = "lightgray"
        elif target_i in ROC_COLORS:
            color = ROC_COLORS[target_i]
        else:
            assert False, f"Target {target_i} not in ROC_COLORS"

        legend_label = targets_df["identifier"][ti]

        # rename CFH to PEC
        if legend_label == "CFH":
            legend_label = "PEC"

        if target_i == target:
            # plot main line thicker
            plt.plot(recall, precision, color=color,
                     lw=3, label="{} AUPRC={:.3f}".format(legend_label, auprc))
        else:
            plt.plot(recall, precision, color=color,
                     lw=1, label="{} AUPRC={:.3f}".format(legend_label, auprc))

    plt.axhline(y=num_pos / num_total, color="gray", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("PRC, classifying positive vs. negative set")

    # Reorder by AUPRC
    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.argsort(auprcs)[::-1]
    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])

    # Save plot
    plt.savefig(f"{out_dir}/prc_class_cross.pdf", dpi=300)
    plt.close("all")
    metrics_df.loc[f"{target}_auprc_class"] = auprcs

    return metrics_df


def get_imb_from_sad(sad_h5: h5py.File, ti: int) -> np.ndarray:
    """
    Calculate allelic imbalance from SAD file and target index.
    """
    ref = sad_h5["REF"][:, :, ti].squeeze()
    alt = sad_h5["ALT"][:, :, ti].squeeze()
    return ref / (ref + alt)

if __name__ == "__main__":
    main()
