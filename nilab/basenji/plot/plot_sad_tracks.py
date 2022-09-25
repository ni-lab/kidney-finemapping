#!/usr/bin/env python

import os
from optparse import OptionParser
from typing import List

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typeguard import typechecked

from nilab.kidney_utils import PLOT_KIDNEY_TARGET_INDICES, KIDNEY_CMAP


def main():
    usage = 'usage: %prog [options] <sad_file>'
    parser = OptionParser(usage)
    parser.add_option('--stats', dest='sad_stats', default='SAD,REF,ALT', help='Comma-separated list of stats to plot')
    parser.add_option('-t', dest='targets_file', default=None, type='str', help='Targets file')
    parser.add_option('-o', dest='out_dir', default='sad_pos_shifts_plots', help='Output directory for script')
    options, args = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    sad_file = args[0]
    os.makedirs(options.out_dir, exist_ok=True)
    options.sad_stats = options.sad_stats.split(',')

    # Read in files
    sad_h5 = h5py.File(sad_file, mode='r')
    targets = pd.read_table(options.targets_file, index_col=0)

    #################################################################
    # Calculate metrics and sort by max sad ratio
    num_snps = sad_h5['snp'].shape[0]
    metrics = ['max_sad_ratio', 'max_sad', 'max_sad_ratio_per_target', 'max_sad_per_target']
    metrics_df = []
    for i in range(num_snps):
        snp_metrics = calc_snp_metrics(i, sad_h5, targets, metrics)
        metrics_df.append(snp_metrics)
    metrics_df = pd.concat(metrics_df, axis=0)
    metrics_df_sorted = metrics_df.sort_values(by="max_sad_ratio", ascending=False)
    metrics_df_sorted.to_csv(f"{options.out_dir}/metrics.csv", sep='\t', index=False, header=True)

    # Plot SAD tracks and save to out directory
    plot_snp_sad_comparison(sad_h5, targets, options.out_dir, options.sad_stats)


@typechecked
def plot_snp_sad_comparison(sad_h5: h5py.File, targets: pd.DataFrame, out_dir: str, sad_stats: List[str]) -> None:
    """
    Plots a comparison of SAD tracks between classification and regression models for a single SNP for each target.
    """

    num_snps = sad_h5['snp'].shape[0]
    num_preds = sad_h5['REF'].shape[0]
    num_pos_per_snp = num_preds // num_snps
    mid = num_pos_per_snp // 2  # seq_len - seq_len // 2 = seq_len // 2 for even seq_len

    num_targets = targets.shape[0]

    for i in range(num_snps):
        ref_allele = sad_h5["ref_allele"][i].decode("utf-8")
        alt_allele = sad_h5["alt_allele"][i].decode("utf-8")
        rsid = sad_h5['snp'][i].decode('utf-8')
        if len(ref_allele) > 1 or len(alt_allele) > 1:
            # skip indels for plotting
            print(f"Skipping indel for plotting: {rsid}")
            continue

        fig, axs = plt.subplots(num_targets, len(sad_stats), figsize=(50, 100))
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
                ax.set_ylabel('{} Score'.format(sad_stat), fontsize=50, labelpad=20)

                # Make ylim scale to 1.1*max y_value across all targets.
                if sad_stat == 'SAD':
                    ax.set_ylim(-y_max * 1.1 / 2, y_max * 1.1 / 2)
                else:
                    ax.set_ylim(y_max * -0.05, y_max * 1.05)
                ax.tick_params(axis="y", direction="inout", length=18, width=3, color="black", pad=15, labelsize=30)
                ax.set_title('{}, {}'.format(rsid, targets['identifier'].iloc[ti]),
                             fontsize=50, pad=15)
                ax.plot(xs, track_ti, color=KIDNEY_CMAP(ti))
                ax.fill_between(xs, track_ti, color=KIDNEY_CMAP(ti))

        fig.tight_layout(pad=4.0)
        plt.savefig(os.path.join(out_dir, '{}.pdf'.format(rsid)), dpi=400)
        plt.close('all')


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

    num_snps = sad_h5['snp'].shape[0]
    num_preds = sad_h5['REF'].shape[0]
    num_pos_per_snp = num_preds // num_snps
    start, end = num_pos_per_snp * i, num_pos_per_snp * (i + 1)

    metric_names = []

    for metric in metrics:
        if metric == 'max_sad':
            snp_metric = np.max(np.abs(sad_h5['SAD'][start:end]))
            metric_names.append(metric)
            snp_metrics.append(snp_metric)
        elif metric == 'max_sad_ratio':
            snp_metric = np.max(np.max(np.abs(sad_h5['SAD'][start:end]), axis=0) / np.max((sad_h5['REF'][start:end]
                                                                                           + sad_h5['ALT'][
                                                                                             start:end]) / 2,
                                                                                          axis=(0, 1)))
            metric_names.append(metric)
            snp_metrics.append(snp_metric)

        elif metric == 'max_sad_per_target':
            snp_metric = np.max(np.abs(sad_h5['SAD'][start:end]), axis=0)

            for ti in PLOT_KIDNEY_TARGET_INDICES:
                metric_names.append('max_sad_{}'.format(targets['identifier'].iloc[ti]))
                snp_metrics.append(snp_metric[ti])

        elif metric == 'max_sad_ratio_per_target':
            snp_metric = np.max(np.abs(sad_h5['SAD'][start:end]), axis=0) / np.max((sad_h5['REF'][start:end] +
                                                                                    sad_h5['ALT'][start:end]) / 2,
                                                                                   axis=(0, 1))
            for ti in PLOT_KIDNEY_TARGET_INDICES:
                metric_names.append('max_sad_ratio_{}'.format(targets['identifier'].iloc[ti]))
                snp_metrics.append(snp_metric[ti])
        else:
            assert False, 'Unrecognized metric {}'.format(metric)

    snp_metrics = pd.DataFrame([[sad_h5['snp'][i].decode('utf-8')] + snp_metrics], columns=['rsid', *metric_names])

    return snp_metrics


if __name__ == '__main__':
    main()
