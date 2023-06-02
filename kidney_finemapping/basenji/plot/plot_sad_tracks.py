#!/usr/bin/env python

import os
from optparse import OptionParser
from typing import List, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typeguard import typechecked
from scipy.interpolate import make_interp_spline

from kidney_finemapping.kidney_utils import PLOT_KIDNEY_TARGET_INDICES, KIDNEY_CMAP


def main():
    """
    Plot SAD tracks for each SNP from the output of compute_sad_shifts.py and save to out directory.
    Shifts reflect prediction window position relative to the SNP, so a shift of -1 means we are predicting one position to the left of the SNP.
    - args[0] <sad_file> - HDF5 file containing SAD shift scores for each SNP
    """
    usage = "usage: %prog [options] <sad_file>"
    parser = OptionParser(usage)
    parser.add_option("-t", dest="targets_file", default=None, type="str", help="Targets file")
    parser.add_option("--overlay", dest="overlay", default=False, action="store_true", help="Overlay ref and alt for plotting")
    parser.add_option("--overlay_lines_only", dest="overlay_lines_only", default=False, action="store_true",
                      help="If using overlay, plot lines only (do not fill between).")
    parser.add_option("-o", dest="out_dir", default="sad_pos_shifts_plots", help="Output directory for script")
    options, args = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    sad_file = args[0]
    os.makedirs(options.out_dir, exist_ok=True)

    # Read in files
    sad_h5 = h5py.File(sad_file, mode="r")
    targets = pd.read_table(options.targets_file, index_col=0)

    #################################################################
    # Calculate metrics and sort by max sad ratio
    num_snps = sad_h5["snp"].shape[0]
    metrics = ["max_sad_ratio", "max_sad", "max_sad_ratio_per_target", "max_sad_per_target"]
    metrics_df = []
    for i in range(num_snps):
        snp_metrics = calc_snp_metrics(i, sad_h5, targets, metrics)
        metrics_df.append(snp_metrics)
    metrics_df = pd.concat(metrics_df, axis=0)
    metrics_df_sorted = metrics_df.sort_values(by="max_sad_ratio", ascending=False)
    metrics_df_sorted.to_csv(f"{options.out_dir}/metrics.csv", sep="\t", index=False, header=True)

    # Plot SAD tracks and save to out directory
    if options.overlay:
        plot_snp_sad_comparison_overlay(sad_h5, targets, options.out_dir, lines_only=options.overlay_lines_only)
    else:
        plot_snp_sad_comparison(sad_h5, targets, options.out_dir, sad_stats=["SAD", "REF", "ALT"])


@typechecked
def plot_snp_sad_comparison(sad_h5: h5py.File, targets: pd.DataFrame, out_dir: str, sad_stats: List[str]) -> None:
    """
    Plots SAD, REF, and ALT tracks for all snps.
    - If overlay_ref_alt is True, we overlay the REF and ALT tracks, using a darker shade for the less accessible allele
    and a lighter shade for the more accessible allele.
    """
    num_snps = sad_h5["snp"].shape[0]
    num_preds = sad_h5["REF"].shape[0]
    num_pos_per_snp = num_preds // num_snps
    mid = num_pos_per_snp // 2  # seq_len - seq_len // 2 = seq_len // 2 for even seq_len

    num_targets = targets.shape[0]

    for i in range(num_snps):
        ref_allele = sad_h5["ref_allele"][i].decode("utf-8")
        alt_allele = sad_h5["alt_allele"][i].decode("utf-8")
        rsid = sad_h5["snp"][i].decode("utf-8")
        if len(ref_allele) > 1 or len(alt_allele) > 1:
            # skip indels for plotting
            print(f"Skipping indel for plotting: {rsid}")
            continue

        # figure dimensions
        width, height = 50, 100
        fig, axs = plt.subplots(num_targets, len(sad_stats), figsize=(width, height))

        xs = np.linspace(-mid, num_pos_per_snp - mid, num_pos_per_snp, endpoint=False)
        xtick_labels = np.linspace(-mid, num_pos_per_snp - mid, num=7, endpoint=True)

        # Calculate ylim across all plots
        y_max = np.max(np.abs(
            [sad_h5[sad_stat][num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze() for sad_stat in sad_stats]))

        for si, sad_stat in enumerate(sad_stats):
            # Flip track to convert SAD pos shifts from prediction-centered to variant-centered tracks
            tracks = np.flip(sad_h5[sad_stat][num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze(), axis=0)

            for plot_i, ti in enumerate(PLOT_KIDNEY_TARGET_INDICES):
                track_ti = tracks[:, ti]
                ax = axs[plot_i, si]

                ax.set_xticks(xtick_labels)
                ax.tick_params(axis="x", direction="inout", length=18, width=3, color="black", pad=15, labelsize=30)
                ax.set_ylabel("{} Score".format(sad_stat), fontsize=50, labelpad=20)

                # Make ylim scale to 1.1*max y_value across all targets.
                if sad_stat == "SAD":
                    ax.set_ylim(-y_max * 1.1 / 2, y_max * 1.1 / 2)
                else:
                    ax.set_ylim(y_max * -0.05, y_max * 1.05)
                ax.tick_params(axis="y", direction="inout", length=18, width=3, color="black", pad=15, labelsize=30)
                ax.set_title("{}, {}".format(rsid, targets["identifier"].iloc[ti]),
                             fontsize=50, pad=15)
                ax.plot(xs, track_ti, color=KIDNEY_CMAP(ti))
                ax.fill_between(xs, track_ti, color=KIDNEY_CMAP(ti))

        fig.tight_layout(pad=4.0)
        plt.savefig(os.path.join(out_dir, "{}.pdf".format(rsid)), dpi=400)
        plt.close("all")


@typechecked
def plot_snp_sad_comparison_overlay(sad_h5: h5py.File, targets: pd.DataFrame, out_dir: str,
                                    lines_only=False) -> None:
    """
    Plots SAD, REF, and ALT tracks for all snps, overlaying the REF and ALT tracks.
    We use a darker shade for the less accessible allele and a lighter shade for the more accessible allele.
    """
    num_snps = sad_h5["snp"].shape[0]
    num_preds = sad_h5["REF"].shape[0]
    num_pos_per_snp = num_preds // num_snps
    mid = num_pos_per_snp // 2  # seq_len - seq_len // 2 = seq_len // 2 for even seq_len

    num_targets = targets.shape[0]

    for i in range(num_snps):
        ref_allele = sad_h5["ref_allele"][i].decode("utf-8")
        alt_allele = sad_h5["alt_allele"][i].decode("utf-8")
        rsid = sad_h5["snp"][i].decode("utf-8")

        if len(ref_allele) > 1 or len(alt_allele) > 1:
            # skip indels for plotting
            print(f"Skipping indel for plotting: {rsid}")
            continue

        # figure dimensions
        width, height = 33, 100
        fig, axs = plt.subplots(num_targets, 2, figsize=(width, height))

        xs = np.linspace(-mid, num_pos_per_snp - mid, num_pos_per_snp, endpoint=False)
        xtick_labels = np.linspace(-mid, num_pos_per_snp - mid, num=7, endpoint=True)

        # Calculate ylim across all plots
        y_max = np.max(np.abs(
            [sad_h5[sad_stat][num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze() for sad_stat in ["SAD", "REF", "ALT"]]))

        for plot_i, ti in enumerate(PLOT_KIDNEY_TARGET_INDICES):
            # For overlaying REF and ALT, put more accessible allele first as measured by mean
            allele_order = sorted(["REF", "ALT"], key=lambda x: np.mean(
                sad_h5[x][num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze()
            ), reverse=True)

            sad_stats = ["SAD"] + allele_order  # order to plot SAD stats in so that we plot the more accessible allele first

            for si, sad_stat in enumerate(sad_stats):
                # Flip track to convert SAD pos shifts from prediction-centered to variant-centered tracks
                tracks = np.flip(sad_h5[sad_stat][num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze(), axis=0)

                track_ti = tracks[:, ti]

                if sad_stat == "SAD":
                    ax = axs[plot_i, 0]
                    ax.set_ylabel("SAD Score", fontsize=50, labelpad=20)
                else:
                    # REF or ALT
                    ax = axs[plot_i, 1]
                    ax.set_ylabel(f"REF (darker)\nALT (lighter)",
                                  fontsize=40, labelpad=20)

                ax.set_xticks(xtick_labels)
                ax.tick_params(axis="x", direction="inout", length=18, width=3, color="black", pad=15, labelsize=30)

                # Make ylim scale to 1.1*max y_value across all targets.
                # Important: must be consistent between SAD tracks and REF / ALT for scaling reasons
                ylim_scale = 1.1
                if sad_stat == "SAD":
                    ax.set_ylim(-y_max * ylim_scale / 2, y_max * ylim_scale / 2)
                else:
                    s = (ylim_scale - 1) / 2
                    ax.set_ylim(y_max * -s, y_max * (1 + s))
                ax.tick_params(axis="y", direction="inout", length=18, width=3, color="black", pad=15, labelsize=30)
                ax.set_title("{}, {}".format(rsid, targets["identifier"].iloc[ti]),
                             fontsize=50, pad=15)

                r, g, b, a = KIDNEY_CMAP(ti)
                if sad_stat == "ALT":
                    # plot alternate allele with lighter color
                    # tinting: https://stackoverflow.com/questions/6615002/given-an-rgb-value-how-do-i-create-a-tint-or-shade
                    if lines_only:
                        # tint much lighter if only plotting the lines
                        tint_factor = 0.7
                    else:
                        tint_factor = 0.5
                    rt = r + (tint_factor * (1 - r))
                    gt = g + (tint_factor * (1 - g))
                    bt = b + (tint_factor * (1 - b))
                    color = (rt, gt, bt, a)
                else:
                    color = (r, g, b, a)

                xnew, track_smooth = smooth_track(xs, track_ti, k=3)

                if lines_only:
                    ax.plot(xnew, track_smooth, color=color, lw=5)
                else:
                    if si == 1:
                        # plot only the top of the more accessible allele
                        less_accessible_sad_stat = sad_stats[2]
                        track_less_accessible = np.flip(sad_h5[less_accessible_sad_stat]
                                                        [num_pos_per_snp * i:num_pos_per_snp * (i + 1)].squeeze(), axis=0)[:, ti]

                        _, track_less_accessible_smooth = smooth_track(xs, track_less_accessible, k=3)
                        ax.fill_between(xnew, track_smooth, track_less_accessible_smooth, color=color)
                    else:
                        # plot SAD and less accessible allele
                        ax.fill_between(xnew, track_smooth, color=color)

        fig.tight_layout(pad=4.0)
        plt.savefig(os.path.join(out_dir, "{}.pdf".format(rsid)), dpi=400)
        plt.close("all")



@typechecked
def calc_snp_metrics(i: int, sad_h5: h5py.File, targets: pd.DataFrame, metrics: List[str]):
    """
    Calculates metrics on SAD tracks for each SNP.

    Args:
     - i: index of SNP in H5 file.
     - sad_h5: SAD scores from basenji_sad_pos_shifts output
     - targets: Dataframe of targets
     - metrics: metrics to calculate
         - max_sad: Max absolute SAD score across track over all targets.
         - max_sad_ratio_per_target: Max (absolute SAD / max ((REF + ALT)/2) across track) for each target.
         - max_sad_ratio: Max over all targets of max_sad_ratio_per_target

    Output:
     - snp_metrics: list with corresponding metrics computed for given SNP
    """
    snp_metrics = []

    num_snps = sad_h5["snp"].shape[0]
    num_preds = sad_h5["REF"].shape[0]
    num_pos_per_snp = num_preds // num_snps
    start, end = num_pos_per_snp * i, num_pos_per_snp * (i + 1)

    metric_names = []

    for metric in metrics:
        if metric == "max_sad":
            snp_metric = np.max(np.abs(sad_h5["SAD"][start:end]))
            metric_names.append(metric)
            snp_metrics.append(snp_metric)
        elif metric == "max_sad_ratio":
            snp_metric = np.max(np.max(np.abs(sad_h5["SAD"][start:end]), axis=0) / np.max((sad_h5["REF"][start:end]
                                                                                           + sad_h5["ALT"][
                                                                                             start:end]) / 2,
                                                                                          axis=(0, 1)))
            metric_names.append(metric)
            snp_metrics.append(snp_metric)

        elif metric == "max_sad_per_target":
            snp_metric = np.max(np.abs(sad_h5["SAD"][start:end]), axis=0)

            for ti in PLOT_KIDNEY_TARGET_INDICES:
                metric_names.append("max_sad_{}".format(targets["identifier"].iloc[ti]))
                snp_metrics.append(snp_metric[ti])

        elif metric == "max_sad_ratio_per_target":
            snp_metric = np.max(np.abs(sad_h5["SAD"][start:end]), axis=0) / np.max((sad_h5["REF"][start:end] +
                                                                                    sad_h5["ALT"][start:end]) / 2,
                                                                                   axis=(0, 1))
            for ti in PLOT_KIDNEY_TARGET_INDICES:
                metric_names.append("max_sad_ratio_{}".format(targets["identifier"].iloc[ti]))
                snp_metrics.append(snp_metric[ti])
        else:
            assert False, "Unrecognized metric {}".format(metric)

    snp_metrics = pd.DataFrame([[sad_h5["snp"][i].decode("utf-8")] + snp_metrics], columns=["rsid", *metric_names])

    return snp_metrics


def smooth_track(xs: np.array, track: np.array, k: int = 3, num_points: int = 50) -> Tuple[np.array, np.array]:
    """
    Smooth the provided track.
    Args:
        - xs: x-values
        - track: y-values
        - k: B-spline degree
        - num_points: number of points to smooth with

    Output:
        - tuple of (smoothed x values, smoothed y values)
    """
    spl = make_interp_spline(xs, track, k=k)
    xnew = np.linspace(xs.min(), xs.max(), num_points)
    track_smooth = spl(xnew)
    return xnew, track_smooth


if __name__ == "__main__":
    main()
