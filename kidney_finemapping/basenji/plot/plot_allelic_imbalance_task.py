# !/usr/bin/env python
from optparse import OptionParser
import h5py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import bisect
from collections import defaultdict

from roc_comparison.compare_auc_delong_xu import delong_roc_variance
from basenji.kidney_utils import NON_TUBULE_EPITHELIAL, ROC_COLORS


'''
allelic_imbalance_plot.py

Given allelic imbalance variant sets and SAD files, output ROC curves and other metrics.
'''


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <variant_set_dir>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir',
                      default='plots'),
    parser.add_option('--neg_mult', dest='neg_mult',
                      default=1, type='int',
                      help='Multiplier for negative set size relative to positive set [Default: %default]')
    parser.add_option('--thresh', dest='thresh',
                      default=0.05, type='float',
                      help='Threshold used for creating allelic imbalance sets')
    parser.add_option('--targets', dest='targets',
                      default='LOH,PT',
                      help='Comma-separated list denoting targets to include in analysis [Default: %default]')
    parser.add_option('-t', dest='targets_file',
                      default=None, type='str',
                      help='File specifying target indexes and labels in table format')
    parser.add_option('--plot_combined', dest='plot_combined',
                      default=False, action='store_true',
                      help='Plot combined reads using PT/LOH predictions.')

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error('Must provide variant sets directory')
    else:
        variant_sets_dir = args[0]

    os.makedirs(options.out_dir, exist_ok=True)
    targets = options.targets.split(',')
    neg_mult = options.neg_mult
    thresh = options.thresh


    #######################################################
    # Load files
    pos_dirs = [os.path.join(variant_sets_dir, '{}_neg{}x_q{}/pos_sad/all_chrs/'.format(target, neg_mult, thresh))
                for target in targets]
    neg_dirs = [os.path.join(variant_sets_dir, '{}_neg{}x_q{}/neg_sad/all_chrs/'.format(target, neg_mult, thresh))
                for target in targets]

    pos_vcf_files = [os.path.join(dir_, 'all_chrs_snps.vcf') for dir_ in pos_dirs]
    neg_vcf_files = [os.path.join(dir_, 'all_chrs_snps.vcf') for dir_ in neg_dirs]

    pos_sad_files = [os.path.join(dir_, 'sad.h5') for dir_ in pos_dirs]
    neg_sad_files = [os.path.join(dir_, 'sad.h5') for dir_ in neg_dirs]

    col_names = names = ['chr', 'pos', 'index', 'ref', 'alt', 'ref_reads', 'alt_reads', 'Pvalue', 'imbalance',
                         'total_reads', 'Qvalue']
    pos_vcfs = [pd.read_table(vcf, sep='\t', names=col_names, index_col=None) for vcf in pos_vcf_files]
    neg_vcfs = [pd.read_table(vcf, sep='\t', names=col_names, index_col=None) for vcf in neg_vcf_files]

    pos_sads = [h5py.File(sad, mode='r') for sad in pos_sad_files]
    neg_sads = [h5py.File(sad, mode='r') for sad in neg_sad_files]

    data_per_target = list(zip(targets, pos_vcfs, neg_vcfs, pos_sads, neg_sads))

    targets_df = pd.read_csv(options.targets_file, sep='\t', index_col=0)

    #######################################################
    # Plot read count distributions for positive and negative sets
    for target, pos_vcf, neg_vcf in zip(targets, pos_vcfs, neg_vcfs):
        plt.title('Read counts for pos and neg sets, {}'.format(target))
        x_max = np.max(pos_vcf['total_reads']) + 1
        plt.hist(pos_vcf['total_reads'], alpha=0.5, bins=np.arange(0, x_max, (x_max // 20) + 1), density=True,
                 label='Positive')
        plt.hist(neg_vcf['total_reads'], alpha=0.5, bins=np.arange(0, x_max, (x_max // 20) + 1), density=True,
                 label='Negative')
        plt.legend()
        plt.savefig(os.path.join(options.out_dir, '{}_count_hist.pdf'.format(target)), dpi=600)
        plt.show()
        plt.close('all')

    #######################################################
    # Main plots: plot ROCs for direction and classification
    # Plot params
    target_indices = [2, 4, 8, 9, 7, 3, 6, 1, 5, 0]
    cmap = plt.get_cmap('tab10')
    figsize = (5.5, 5)
    baseline_color = 'gray'
    lw_main = 2

    # Plot AI direction
    plt.figure(figsize=figsize)
    metrics = defaultdict(list)
    for data in data_per_target:
        target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
        ti = np.where(targets_df['identifier'] == target)[0].squeeze()
        imb_labels = pos_vcf['imbalance'] > 0.5
        imb_preds = 2 ** pos_sad['IMB'][:, ti]

        fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
        roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
        metrics[target].append(np.sqrt(auc_cov))
        plt.plot(fpr, tpr, color=cmap(ti),
                 lw=lw_main, label='{} AUROC={:.3f}'.format(target, roc_auc))
        plt.legend()

    # Hard coded filepaths
    if options.plot_combined:
        pos_vcf = pd.read_table(
            os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/pos_sad/all_chrs/all_chrs_snps.vcf'),
            sep='\t', names=col_names, index_col=None)
        neg_vcf = pd.read_table(
            os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/neg_sad/all_chrs/all_chrs_snps.vcf'),
            sep='\t', names=col_names, index_col=None)
        pos_sad = h5py.File(os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/pos_sad/all_chrs/sad.h5'), mode='r')
        neg_sad = h5py.File(os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/neg_sad/all_chrs/sad.h5'), mode='r')
        comb_cmap = plt.get_cmap('Set3')

        for target in ['PT', 'LOH', 'DT']:
            ti = np.where(targets_df['identifier'] == target)[0].squeeze()
            imb_labels = pos_vcf['imbalance'] > 0.5
            imb_preds = 2 ** pos_sad['IMB'][:, ti]
            fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
            roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
            metrics['combined_{}'.format(target)].append(np.sqrt(auc_cov))
            plt.plot(fpr, tpr, color=comb_cmap(ti),
                     lw=lw_main, label='combined_{} AUROC={:.3f}'.format(target, roc_auc))
            plt.legend()

    plt.plot([0, 1], [0, 1], color=baseline_color, lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right", framealpha=1)
    plt.title('Predicting direction of allelic imbalance')
    plt.savefig(os.path.join(options.out_dir, 'all_roc_direc.pdf'), dpi=600)
    plt.show()

    # Classification
    target_indices = [2, 4, 8, 9, 7, 3, 6, 1, 5, 0]
    cmap = plt.get_cmap('tab10')

    plt.figure(figsize=figsize)
    metrics = defaultdict(list)
    for data in data_per_target:
        target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
        ti = np.where(targets_df['identifier'] == target)[0].squeeze()
        imb_labels = np.concatenate([np.ones(pos_vcf['imbalance'].shape[0]),
                                     np.zeros(neg_vcf['imbalance'].shape[0])])
        pos_imbs = 2 ** pos_sad['IMB'][:, ti]
        neg_imbs = 2 ** neg_sad['IMB'][:, ti]
        imb_preds = np.concatenate([pos_imbs, neg_imbs])
        imb_preds = np.abs(imb_preds - 0.5)  # Classify symmetrically about 0.5

        fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)  # Positive vs. negative
        roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
        metrics[target].append(np.sqrt(auc_cov))
        plt.plot(fpr, tpr, color=cmap(ti),
                 lw=lw_main, label='{} AUROC={:.3f}'.format(target, roc_auc))
        plt.legend()

    if options.plot_combined:
        pos_vcf = pd.read_table(
            os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/pos_sad/all_chrs/all_chrs_snps.vcf'),
            sep='\t', names=col_names, index_col=None)
        neg_vcf = pd.read_table(
            os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/neg_sad/all_chrs/all_chrs_snps.vcf'),
            sep='\t', names=col_names, index_col=None)
        pos_sad = h5py.File(os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/pos_sad/all_chrs/sad.h5'), mode='r')
        neg_sad = h5py.File(os.path.join(variant_sets_dir, 'combined_neg2x_q0.01/neg_sad/all_chrs/sad.h5'), mode='r')
        comb_cmap = plt.get_cmap('Set3')

        for target in ['PT', 'LOH', 'DT']:
            ti = np.where(targets_df['identifier'] == target)[0].squeeze()
            imb_labels = np.concatenate([np.ones(pos_vcf['imbalance'].shape[0]),
                                         np.zeros(neg_vcf['imbalance'].shape[0])])
            pos_imbs = 2 ** pos_sad['IMB'][:, ti]
            neg_imbs = 2 ** neg_sad['IMB'][:, ti]
            imb_preds = np.concatenate([pos_imbs, neg_imbs])
            imb_preds = np.abs(imb_preds - 0.5)  # Classify symmetrically about 0.5
            fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
            roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
            metrics['combined_{}'.format(target)].append(np.sqrt(auc_cov))
            plt.plot(fpr, tpr, color=comb_cmap(ti),
                     lw=lw_main, label='combined_{} AUROC={:.3f}'.format(target, roc_auc))

    plt.plot([0, 1], [0, 1], color=baseline_color, lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right", framealpha=1)
    plt.title('Predicting allelic imbalance')
    plt.savefig(os.path.join(options.out_dir, 'all_roc_cls.pdf'), dpi=600)
    plt.show()

    #######################################################
    # Plot ROCs for cross cell type predictions
    # Predict direction of allelic imbalance in positive set
    metrics_df = pd.DataFrame(columns=targets_df['identifier'][target_indices])
    for data in data_per_target:
        target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
        imb_labels = pos_vcf['imbalance'] > 0.5

        auroc_direc = []  # keep track of auroc to sort legend
        std_direc = []
        fpr_thresh = 0.05

        plt.figure(figsize=figsize)

        for ti in target_indices:
            target_i = targets_df['identifier'][ti]
            imb_preds = 2 ** pos_sad['IMB'][:, ti]
            fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
            roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
            auroc_direc.append(roc_auc)
            std_direc.append(np.sqrt(auc_cov))

            if target_i in NON_TUBULE_EPITHELIAL:
                color = "lightgray"
            elif target_i in ROC_COLORS:
                color = ROC_COLORS[target_i]
            else:
                color = cmap(ti)

            legend_label = targets_df['identifier'][ti]

            # rename CFH to PEC
            if legend_label == "CFH":
                legend_label = "PEC"

            if target_i == target:
                plt.plot(fpr, tpr, color=color,
                         lw=3, label='{} AUROC={:.3f}'.format(legend_label, roc_auc))
                th_i = bisect.bisect(fpr, fpr_thresh) - 1
                fpr_th = fpr[th_i]
                tpr_th = tpr[th_i]
                ai_th = thresholds[th_i]
            else:
                plt.plot(fpr, tpr, color=color,
                         lw=1, label='{} AUROC={:.3f}'.format(legend_label, roc_auc))
        plt.plot([0, 1], [0, 1], color=baseline_color, lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC, predicting allelic imbalance direction \n (FPR: {:.3f}, TPR: {:.3f}, AI: {:.3f})'.format(fpr_th,
                                                                                                                 tpr_th,
                                                                                                                 ai_th))
        plt.legend(loc="lower right", framealpha=1)
        plt.legend(loc="lower right", framealpha=1)

        # reorder by auc
        handles, labels = plt.gca().get_legend_handles_labels()
        order = np.argsort(auroc_direc)[::-1]  # sort in descending order of auroc
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])

        plt.savefig(os.path.join(options.out_dir, '{}_roc_direc_cross.pdf'.format(target)), dpi=600)
        plt.show()
        plt.close('all')
        metrics_df.loc['{}_auroc_direc'.format(target)] = auroc_direc
        metrics_df.loc['{}_std_direc'.format(target)] = std_direc

    #######################################################
    # Discriminating positive set from negative set
    for data in data_per_target:
        target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
        imb_labels = np.concatenate([np.ones(pos_vcf['imbalance'].shape[0]),
                                     np.zeros(neg_vcf['imbalance'].shape[0])])
        auroc_cls = []
        std_cls = []

        plt.figure(figsize=figsize)

        for ti in target_indices:
            target_i = targets_df['identifier'][ti]
            pos_imbs = 2 ** pos_sad['IMB'][:, ti]
            neg_imbs = 2 ** neg_sad['IMB'][:, ti]
            imb_preds = np.concatenate([pos_imbs, neg_imbs])
            imb_preds = np.abs(imb_preds - 0.5)  # Classify symmetrically about 0.5
            fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)  # Positive vs. negative
            roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
            auroc_cls.append(roc_auc)
            std_cls.append(np.sqrt(auc_cov))

            if target_i in NON_TUBULE_EPITHELIAL:
                color = "lightgray"
            elif target_i in ROC_COLORS:
                color = ROC_COLORS[target_i]
            else:
                color = cmap(ti)

            legend_label = targets_df['identifier'][ti]

            # rename CFH to PEC
            if legend_label == "CFH":
                legend_label = "PEC"

            if target == target_i:
                plt.plot(fpr, tpr, color=color,
                         lw=3, label='{} AUROC={:.3f}'.format(legend_label, roc_auc))
                th_i = bisect.bisect(fpr, fpr_thresh) - 1
                fpr_th = fpr[th_i]
                tpr_th = tpr[th_i]
                ai_th = thresholds[th_i]
                plt.hlines(tpr_th, 0, fpr_th, linestyle='dotted', color='orange', lw=1, zorder=100)
                plt.vlines(fpr_th, 0, tpr_th, linestyle='dotted', color='orange', lw=1, zorder=100)
            else:
                plt.plot(fpr, tpr, color=color,
                         lw=1, label='{} AUROC={:.3f}'.format(legend_label, roc_auc))

        plt.plot([0, 1], [0, 1], color=baseline_color, lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(
            'ROC, classifying positive vs. negative set \n (FPR: {:.3f}, TPR: {:.3f}, AI neg interval: ({:.3f}, {:.3f})'.format(
                fpr_th, tpr_th, 0.5 - ai_th, 0.5 + ai_th))

        plt.legend(loc="lower right", framealpha=1)

        # reorder by auc
        handles, labels = plt.gca().get_legend_handles_labels()
        order = np.argsort(auroc_cls)[::-1]  # sort in descending order of auroc
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])

        plt.savefig(os.path.join(options.out_dir, '{}_roc_cls_cross.pdf'.format(target)), dpi=600)
        plt.show()
        plt.close('all')

        metrics_df.loc['{}_auroc_cls'.format(target)] = auroc_cls
        metrics_df.loc['{}_std_cls'.format(target)] = std_cls

        # #######################################################
        # # AUPRC
        # #######################################################
        # Same cell type AUPRCs
        #######################################################
        # Plot ROCs
        # Classification
        target_indices = [2, 4, 8, 9, 7, 3, 6, 1, 5, 0]
        cmap = plt.get_cmap('tab10')

        plt.figure(figsize=figsize)
        metrics = defaultdict(list)
        for data in data_per_target:
            target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
            ti = np.where(targets_df['identifier'] == target)[0].squeeze()
            imb_labels = np.concatenate([np.ones(pos_vcf['imbalance'].shape[0]),
                                         np.zeros(neg_vcf['imbalance'].shape[0])])
            pos_imbs = 2 ** pos_sad['IMB'][:, ti]
            neg_imbs = 2 ** neg_sad['IMB'][:, ti]
            imb_preds = np.concatenate([pos_imbs, neg_imbs])
            imb_preds = np.abs(imb_preds - 0.5)  # Classify symmetrically about 0.5
            precision, recall, thresholds = precision_recall_curve(imb_labels, imb_preds)  # Positive vs. negative
            auprc = auc(recall, precision)

            plt.plot(recall, precision, color=cmap(ti),
                     lw=2, label='{} AUPRC={:.3f}'.format(targets_df['identifier'][ti], auprc))
        num_pos = pos_vcf.shape[0]
        num_total = num_pos + neg_vcf.shape[0]
        plt.axhline(y=num_pos / num_total, color=baseline_color, lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(loc="upper right", framealpha=1)
        plt.title('Predicting allelic imbalance')
        plt.savefig(os.path.join(options.out_dir, 'all_prc_class.pdf'), dpi=600)

        #######################################################
        # Cross-cell type AUPRCs
        for data in data_per_target:
            target, pos_vcf, neg_vcf, pos_sad, neg_sad = data
            imb_labels = np.concatenate(
                [np.ones(pos_vcf['imbalance'].shape[0]), np.zeros(neg_vcf['imbalance'].shape[0])])

            auprc_cls = []

            plt.figure(figsize=figsize)

            for ti in target_indices:
                target_i = targets_df['identifier'][ti]
                pos_imbs = 2 ** pos_sad['IMB'][:, ti]
                neg_imbs = 2 ** neg_sad['IMB'][:, ti]
                imb_preds = np.concatenate([pos_imbs, neg_imbs])
                imb_preds = np.abs(imb_preds - 0.5)  # Classify symmetrically about 0.5
                precision, recall, thresholds = precision_recall_curve(imb_labels, imb_preds)  # Positive vs. negative
                auprc = auc(recall, precision)
                auprc_cls.append(auprc)

                if target_i in NON_TUBULE_EPITHELIAL:
                    color = "lightgray"
                elif target_i in ROC_COLORS:
                    color = ROC_COLORS[target_i]
                else:
                    color = cmap(ti)

                legend_label = targets_df['identifier'][ti]

                # rename CFH to PEC
                if legend_label == "CFH":
                    legend_label = "PEC"

                if target == target_i:
                    plt.plot(recall, precision, color=color,
                             lw=3, label='{} AUPRC={:.3f}'.format(legend_label, auprc))
                else:
                    plt.plot(recall, precision, color=color,
                             lw=1, label='{} AUPRC={:.3f}'.format(legend_label, auprc))

            plt.axhline(y=num_pos / num_total, color=baseline_color, lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title('PRC, classifying positive vs. negative set')
            plt.legend(loc="lower left", framealpha=1, prop={'size': 8}, bbox_to_anchor=(1.1, 0))

            # reorder by auc
            handles, labels = plt.gca().get_legend_handles_labels()
            order = np.argsort(auprc_cls)[::-1]  # sort in descending order of auroc
            plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])

            plt.savefig(os.path.join(options.out_dir, '{}_prc_cls_cross.pdf'.format(target)), dpi=600,
                        bbox_inches='tight')
            plt.close('all')

            metrics_df.loc['auprc_cls'] = auprc_cls

    #######################################################
    # Save standard error and AUC metrics to file
    metrics_df.to_csv(os.path.join(options.out_dir, 'metrics.tsv'), sep='\t', index=True, header=True)


if __name__ == '__main__':
    main()
