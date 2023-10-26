#!/bin/bash
# Commands used to generate PoPS scores from kidney function GWAS summary statistics

./magma --annotate \
    --snp-loc /clusterfs/nilah/rkchung/data/gene_prio/g1000_eur.bim \
    --gene-loc /clusterfs/nilah/rkchung/data/gene_prio/gene_loc.txt \
    --out /clusterfs/nilah/rkchung/data/gene_prio/magma_annot

./magma --bfile /clusterfs/nilah/rkchung/data/gene_prio/g1000_eur \
    --gene-annot /clusterfs/nilah/rkchung/data/gene_prio/magma_annot.genes.annot \
    --pval /clusterfs/nilah/rkchung/data/gene_prio/hg19_egfr_creat_cys_DUPSREMOVED.rsid.txt ncol=N \
    --gene-model snp-wise=mean \
    --out /clusterfs/nilah/rkchung/data/gene_prio/magma_hg19_egfr_creat_cys

# PoPS preprocessed features can be downloaded from https://www.dropbox.com/sh/o6t5jprvxb8b500/AADZ8qD6Rpz4uvCk0b5nUnPaa/data or linked from their website https://www.finucanelab.org/data                                                                                                                                                                         
python ~/pops/munge_feature_directory.py  \
    --gene_annot_path /clusterfs/nilah/rkchung/data/gene_prio/gene_loc_tss.txt \
    --feature_dir /clusterfs/nilah/rkchung/data/gene_prio/pops_features/ \
    --save_prefix /clusterfs/nilah/rkchung/data/gene_prio/pops_proc_features  \
    --max_cols 500

python ~/pops/pops.py  --gene_annot_path /clusterfs/nilah/rkchung/data/gene_prio/gene_loc_tss.txt  --feature_mat_prefix /clusterfs/nilah/rkchung/data/gene_prio/pops_proc_features/pops_proc_features  --num_feature_chunks 115  --magma_prefix /clusterfs/nilah/rkchung/data/gene_prio/magma_hg19_egfr_creat_cys --control_features_path /clusterfs/nilah/rkchung/data/gene_prio/control.features --out_prefix /clusterfs/nilah/rkchung/data/gene_prio/PoPS.hg19_egfr_creat_cys.out
