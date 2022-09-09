#!/usr/bin/env python
import os
from optparse import OptionParser

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from nilab.basenji.basenji_utils import plots
from nilab.basenji.basenji_utils.dna_io import dna_1hot


def main():
    usage = 'usage: %prog [options] <ism_summed_pos_shifts_dir> <sad_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gain',
                      default=False, action='store_true',
                      help='Draw a sequence logo for the gain score, too [Default: %default]')
    parser.add_option('-l', dest='mut_len',
                      default=20, type='int',
                      help='Length of center sequence to mutate used in basenji_ism_summed_pos_shifts [Default: %default]')
    parser.add_option('-m', dest='min_limit',
                      default=0.05, type='float',
                      help='Minimum heatmap limit [Default: %default]')
    parser.add_option('-o', dest='out_dir',
                      default='sat_plot', help='Output directory [Default: %default]')
    parser.add_option('--png', dest='save_png',
                      default=False, action='store_true',
                      help='Write PNG as opposed to PDF [Default: %default]')
    parser.add_option('-t', dest='targets_file',
                      default=None, type='str',
                      help='File specifying target indexes and labels in table format')
    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 2
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    ism_summed_pos_shifts_dir = args[0]
    sad_h5_file = args[1]

    os.makedirs(options.out_dir, exist_ok=True)
    save_ext = 'pdf'
    if options.save_png:
        save_ext = 'png'

    # Assume symmetric
    mut_up = mut_down = options.mut_len // 2

    # Set up dataframe with rsids, ref/alt ism seqs, and ref/alt alleles
    ism_info_file = f"{ism_summed_pos_shifts_dir}/info.h5"
    ism_info = h5py.File(ism_info_file, "r")

    sad_h5 = h5py.File(sad_h5_file, "r")
    sad_rsids = [rsid.decode("utf-8") for rsid in sad_h5["snp"]]
    sad_ref_alleles = [allele.decode("utf-8") for allele in sad_h5["ref_allele"]]
    sad_alt_alleles = [allele.decode("utf-8") for allele in sad_h5["alt_allele"]]
    sad_df = pd.DataFrame({"sad_rsid": sad_rsids, "ref_allele": sad_ref_alleles, "alt_allele": sad_alt_alleles})

    rsids = [rsid.decode('utf-8') for rsid in ism_info['rsid']]
    ref_seqs = [ref_seq.decode('utf-8') for ref_seq in ism_info['ref_ism_seq']]
    alt_seqs = [alt_seq.decode('utf-8') for alt_seq in ism_info['alt_ism_seq']]

    snp_seq_df = pd.DataFrame({"rsid": rsids, "ref_seq": ref_seqs, "alt_seq": alt_seqs})
    snp_seq_df = snp_seq_df.merge(sad_df, how="inner", left_on="rsid", right_on="sad_rsid", validate="1:1").drop("sad_rsid", axis=1)
    snp_seq_df = snp_seq_df.set_index("rsid")

    # determine targets
    targets_df = pd.read_table(options.targets_file, index_col=0)
    num_targets = targets_df.shape[0]

    # determine plot region
    alleles = ['ref', 'alt']
    for rsid, row in tqdm(snp_seq_df.iterrows(), total=snp_seq_df.shape[0]):
        for allele in alleles:
            if allele == "ref":
                preds_dir = "ref_ism_preds"
                ism_seq_col = "ref_seq"
                allele_col = "ref_allele"
            else:
                preds_dir = "alt_ism_preds"
                ism_seq_col = "alt_seq"
                allele_col = "alt_allele"

            plot_out_dir = f"{options.out_dir}/{allele}"
            os.makedirs(plot_out_dir, exist_ok=True)
            plot_start = 0
            plot_end = len(row[ism_seq_col])

            # plot attributes
            sns.set(style='white', font_scale=1)
            spp = subplot_params(plot_end - plot_start)

            # For given rsid, compute scores relative to reference/alternate (original sequence before ISM)
            ism_scores = np.load(f"{ism_summed_pos_shifts_dir}/{preds_dir}/{rsid}.npy")
            seq_1hot = dna_1hot(snp_seq_df.loc[rsid, ism_seq_col])
            ref_scores = np.expand_dims(ism_scores[seq_1hot], axis=1)
            scores_rel_ref = ism_scores - ref_scores

            for ti in range(num_targets):
                # compute scores relative to reference
                delta_ti = scores_rel_ref[:, :, ti]

                # compute loss and gain
                delta_loss = delta_ti.min(axis=1)
                delta_gain = delta_ti.max(axis=1)

                # setup plot
                plt.figure(figsize=(20 * (plot_end - plot_start) / 20, 6))  # scale plot based on end and start
                if options.gain:
                    grid_rows = 4
                else:
                    grid_rows = 3
                row_i = 0
                ax_logo_loss = plt.subplot2grid(
                    (grid_rows, spp['heat_cols']), (row_i, spp['logo_start']),
                    colspan=spp['logo_span'])
                row_i += 1
                if options.gain:
                    ax_logo_gain = plt.subplot2grid(
                        (grid_rows, spp['heat_cols']), (row_i, spp['logo_start']),
                        colspan=spp['logo_span'])
                    row_i += 1
                ax_sad = plt.subplot2grid(
                    (grid_rows, spp['heat_cols']), (row_i, spp['sad_start']),
                    colspan=spp['sad_span'])
                row_i += 1
                ax_heat = plt.subplot2grid(
                    (grid_rows, spp['heat_cols']), (row_i, 0), colspan=spp['heat_cols'])

                # plot sequence logo
                plot_seqlogo(ax_logo_loss, seq_1hot, -delta_loss)
                if options.gain:
                    plot_seqlogo(ax_logo_gain, seq_1hot, delta_gain)

                # plot SAD
                plot_sad(ax_sad, delta_loss, delta_gain)

                # plot heat map
                plot_heat(ax_heat, delta_ti.T, options.min_limit, row[allele_col], mut_up)

                plt.tight_layout()
                plt.savefig(f"{plot_out_dir}/{rsid}_{targets_df['identifier'].iloc[ti]}.{save_ext}", dpi=300)
                plt.show()
                plt.close()


def expand_4l(sat_lg_ti, seq_1hot):
    """ Expand

      In:
          sat_lg_ti (l array): Sat mut loss/gain scores for a single sequence and
          target.
          seq_1hot (Lx4 array): One-hot coding for a single sequence.

      Out:
          sat_loss_4l (lx4 array): Score-hot coding?

      """

    # determine satmut length
    satmut_len = sat_lg_ti.shape[0]

    # jump to satmut region in one hot coded sequence
    ssi = int((seq_1hot.shape[0] - satmut_len) // 2)

    # filter sequence for satmut region
    seq_1hot_sm = seq_1hot[ssi:ssi + satmut_len, :]

    # tile loss scores to align
    sat_lg_tile = np.tile(sat_lg_ti, (4, 1)).T

    # element-wise multiple
    sat_lg_4l = np.multiply(seq_1hot_sm, sat_lg_tile)

    return sat_lg_4l


def plot_heat(ax, sat_delta_ti, min_limit, allele, mut_up):
    """ Plot satmut deltas.

      Args:
          ax (Axis): matplotlib axis to plot to.
          sat_delta_ti (4 x L_sm array): Single target delta matrix for saturated mutagenesis region,
          min_limit (float): Minimum heatmap limit.
          allele (string): allele (either ref or alt)
          mut_up: number of nucleotides upstream for ism
      """

    vlim = max(min_limit, np.nanmax(np.abs(sat_delta_ti)))
    yticks = 'ACGT'
    sns.heatmap(
        sat_delta_ti,
        linewidths=0,
        cmap='RdBu_r',
        vmin=-vlim,
        vmax=vlim,
        xticklabels=False,
        ax=ax)

    # Plot stars on variant nucleotides
    x_is = []
    y_is = []
    for i, nt in enumerate(allele):
        x_is.append(mut_up + i + 0.5)
        y_is.append(yticks.find(nt) + 0.5)
    ax.scatter(x_is, y_is, s=125, marker='*', color='black')

    ax.yaxis.set_ticklabels(yticks, rotation='horizontal')


def plot_predictions(ax, preds, satmut_len, seq_len, buffer):
    """ Plot the raw predictions for a sequence and target
          across the specificed saturated mutagenesis region.

      Args:
          ax (Axis): matplotlib axis to plot to.
          preds (L array): Target predictions for one sequence.
          satmut_len (int): Satmut length from which to determine
                             the values to plot.
          seq_len (int): Full sequence length.
          buffer (int): Ignored buffer sequence on each side
      """

    # repeat preds across pool width
    target_pool = (seq_len - 2 * buffer) // preds.shape[0]
    epreds = preds.repeat(target_pool)

    satmut_start = (epreds.shape[0] - satmut_len) // 2
    satmut_end = satmut_start + satmut_len

    ax.plot(epreds[satmut_start:satmut_end], linewidth=1)
    ax.set_xlim(0, satmut_len)
    ax.axhline(0, c='black', linewidth=1, linestyle='--')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)


def plot_sad(ax, sat_loss_ti, sat_gain_ti):
    """ Plot loss and gain SAD scores.

      Args:
          ax (Axis): matplotlib axis to plot to.
          sat_loss_ti (L_sm array): Minimum mutation delta across satmut length.
          sat_gain_ti (L_sm array): Maximum mutation delta across satmut length.
      """

    rdbu = sns.color_palette('RdBu_r', 10)

    ax.plot(-sat_loss_ti, c=rdbu[0], label='loss', linewidth=1)
    ax.plot(sat_gain_ti, c=rdbu[-1], label='gain', linewidth=1)
    ax.set_xlim(0, len(sat_loss_ti))
    ax.legend()
    # ax_sad.grid(True, linestyle=':')

    ax.xaxis.set_ticks([])
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)


def plot_seqlogo(ax, seq_1hot, sat_score_ti, pseudo_pct=0.05):
    """ Plot a sequence logo for the loss/gain scores.

      Args:
          ax (Axis): matplotlib axis to plot to.
          seq_1hot (Lx4 array): One-hot coding of a sequence.
          sat_score_ti (L_sm array): Minimum mutation delta across satmut length.
          pseudo_pct (float): % of the max to add as a pseudocount.
      """
    sat_score_cp = sat_score_ti.copy()
    satmut_len = len(sat_score_ti)

    # add pseudocounts
    sat_score_cp += pseudo_pct * sat_score_cp.max()

    # expand
    sat_score_4l = expand_4l(sat_score_cp, seq_1hot)

    plots.seqlogo(sat_score_4l, ax)


def subplot_params(seq_len):
    """ Specify subplot layout parameters for various sequence lengths. """
    if seq_len < 500:
        spp = {
            'heat_cols': 400,
            'sad_start': 1,
            'sad_span': 321,
            'logo_start': 0,
            'logo_span': 323
        }
    else:
        spp = {
            'heat_cols': 400,
            'sad_start': 1,
            'sad_span': 320,
            'logo_start': 0,
            'logo_span': 322
        }

    return spp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
