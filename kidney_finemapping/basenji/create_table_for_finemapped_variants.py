#!/usr/bin/env python
import os
from optparse import OptionParser
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


################################################################################
# main
################################################################################
def main():
    """
    Given a list of finemapped variants, their SAD scores, and allelic imbalance data, create a table with all relevant
    information for each variant.
    - args[0] <finemap_file> - File containing finemapped variants
    - args[1] <imbalance_data_dir> - Directory containing allelic imbalance data
    - args[2] <sad_file> - File containing SAD scores for finemapped variants
    """
    usage = 'usage: %prog [options] <finemap_file> <imbalance_data_dir> <sad_file>'
    parser = OptionParser(usage)
    parser.add_option('-t', dest='targets_file',
                      default=None, type='str',
                      help='File specifying target indexes and labels in table format')
    parser.add_option('-o', dest='out_dir',
                      default='ukb_table')

    (options, args) = parser.parse_args()
    num_expected_args = 3
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    finemap_file = args[0]
    allelic_imbalance_dir = args[1]
    sad_file = args[2]

    os.makedirs(options.out_dir, exist_ok=True)
    targets = ['CD', 'DT', 'PT', 'LOH', 'Tubule', 'combined']  # Only include scores from these cell types

    targets_df = pd.read_csv(options.targets_file, sep='\t', index_col=0)
    finemap_df = pd.read_table(finemap_file, header=0, index_col=4, sep='\t')

    #######################################################
    # Create table with predicted imbalance for relevant cell types
    sad_h5 = h5py.File(sad_file, mode='r')

    imb_data = {}
    for ti, target in enumerate(targets_df['identifier']):
        if target in targets:
            ref_scores = np.array(sad_h5['REF'][:, :, ti]).squeeze(1)
            alt_scores = np.array(sad_h5['ALT'][:, :, ti]).squeeze(1)
            imb_data['Pred_IMB_{}'.format(target)] = ref_scores / (ref_scores + alt_scores)

    imb_df = pd.DataFrame(imb_data, index=[snp.decode('utf-8') for snp in sad_h5['snp']])
    imb_df['hg38_bp'] = sad_h5['pos']
    imb_df['A1'] = [allele.decode('utf-8') for allele in sad_h5['ref_allele']]
    imb_df['A2'] = [allele.decode('utf-8') for allele in sad_h5['alt_allele']]

    df_imb = finemap_df.merge(imb_df, how='inner', left_index=True, right_index=True, suffixes=('_finemap', ''), validate="1:1")
    cols = df_imb.columns.tolist()
    cols = cols[0:1] + cols[-3:] + cols[4:-3]
    df_imb = df_imb[cols]  # Reorder columns and select ref and alt allele from sad_h5
    df_imb = df_imb.drop_duplicates()  # remove duplicate rows from duplicated rs760077

    #######################################################
    # Add allelic imbalance data

    # convenience column for joining
    df_imb["chr_pos_a1_a2"] = df_imb["hg_38_chrom"] + "_" + \
                              df_imb["hg38_bp"].astype(str) + "_" + \
                              df_imb["A1"] + "_" + \
                              df_imb["A2"]

    df_imb = df_imb.reset_index()

    for target in targets:
        ai_file = f"{allelic_imbalance_dir}/all_{target}q10.tsv"
        ai_df = pd.read_csv(ai_file, sep="\t")
        ai_df = ai_df[['Chr', 'Pos', 'Ref', 'Alt', '#Ref', '#Alt', 'P-value', '#Ref/(#Ref+#Alt)']]  # Pos in hg38
        ai_df["chr_pos_a1_a2"] = "chr" + ai_df["Chr"].astype(str) + "_" + \
                                ai_df["Pos"].astype(str) + "_" + \
                                ai_df["Ref"] + "_" + \
                                ai_df["Alt"]
        ai_df = ai_df.drop(['Chr', 'Pos', 'Ref', 'Alt'], axis=1)
        col_names = (ai_df.columns[:-1] + f"_{target}").tolist() + [ai_df.columns[-1]]  # add target suffix to all except chr_pos_a1_a2
        ai_df.columns = col_names

        df_imb = df_imb.merge(ai_df, how="left", left_on="chr_pos_a1_a2", right_on="chr_pos_a1_a2", validate="1:1")

    df_imb = df_imb.drop(["chr_pos_a1_a2", "A1_finemap", "A2_finemap"], axis=1)
    df_imb = df_imb.set_index("index")

    df_imb.to_csv(Path(f"{options.out_dir}/table.tsv"), sep='\t')


if __name__ == '__main__':
    main()
