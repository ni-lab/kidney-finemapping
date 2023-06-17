#!/usr/bin/env python
import glob
import os
import shutil
from optparse import OptionParser

import h5py
import numpy as np
from natsort import natsorted


def main():
    usage = 'usage: %prog [options] <ism_summed_pos_shifts_dir>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir',
                      default=None)

    (options, args) = parser.parse_args()

    # Setup
    num_expected_args = 1
    if len(args) != num_expected_args:
        parser.error(f"Incorrect number of arguments, expected {num_expected_args} arguments but got {len(args)}")

    ism_dir = args[0]

    os.makedirs(options.out_dir, exist_ok=True)

    # Create h5 containing ISMs from all batches
    ism_files = natsorted(glob.glob(f'{ism_dir}/chr*/info.h5'))
    assert len(ism_files) == 22, 'Expected 22 ism files for 22 chromosomes, got {} files'.format(len(ism_files))

    isms = [h5py.File(ism_file, mode='r') for ism_file in ism_files]

    # merge info files
    info_file_merged = f'{options.out_dir}/info.h5'

    rsids = np.array([rsid for ism in isms for rsid in ism['rsid']], dtype='S')
    with h5py.File(info_file_merged, 'w') as info_h5:
        info_h5.create_dataset('rsid', data=rsids)
        info_h5.create_dataset('ref_ism_seq',
                               data=np.array([ref_ism_seq for ism in isms for ref_ism_seq in ism['ref_ism_seq']], dtype='S'))
        info_h5.create_dataset('alt_ism_seq',
                               data=np.array([alt_ism_seq for ism in isms for alt_ism_seq in ism['alt_ism_seq']], dtype='S'))

    # merge ref ism preds dir
    ref_pred_merged_dir = f'{options.out_dir}/ref_ism_preds'
    os.makedirs(ref_pred_merged_dir, exist_ok=True)
    ref_preds_dirs = natsorted(glob.glob(f'{ism_dir}/chr*/ref_ism_preds'))
    assert len(ref_preds_dirs) == 22, "Expected 22 directories for ref ism preds"

    pred_count = 0
    for ref_pred_dir in ref_preds_dirs:
        for file in glob.glob(f'{ref_pred_dir}/*.npy'):
            shutil.copy(file, ref_pred_merged_dir)
            pred_count += 1

    assert len(rsids) == pred_count, "Number of SNPs does not match number of ref pred files"

    # merge alt ism preds dir
    alt_pred_merged_dir = f'{options.out_dir}/alt_ism_preds'
    os.makedirs(alt_pred_merged_dir, exist_ok=True)
    alt_preds_dirs = natsorted(glob.glob(f'{ism_dir}/chr*/alt_ism_preds'))
    assert len(alt_preds_dirs) == 22, "Expected 22 directories for alt ism preds"

    pred_count = 0
    for alt_pred_dir in alt_preds_dirs:
        for file in glob.glob(f'{alt_pred_dir}/*.npy'):
            shutil.copy(file, alt_pred_merged_dir)
            pred_count += 1

    assert len(rsids) == pred_count, "Number of SNPs does not match number of alt pred files"


if __name__ == '__main__':
    main()
