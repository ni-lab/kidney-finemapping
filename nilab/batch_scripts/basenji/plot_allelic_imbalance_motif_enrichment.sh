conda activate kidney_finemapping

cd /home/rshuai/research/ni-lab/kidney_finemapping/nilab || exit

for TARGET in PT LOH DT; do
  python3 basenji/plot/plot_allelic_imbalance_motif_enrichment.py \
  /home/rshuai/research/ni-lab/analysis/kidney_data/allelic_imbalance/motif_enrichment/tf_pseudobulk_Pseudobulk_Wilson_TF_analysis.csv \
  /home/rshuai/research/ni-lab/analysis/kidney_data/allelic_imbalance/motif_enrichment/${TARGET}_me_variants_neg7x_min_1e-3_diff_1_q0.01_JASPAR2020_CORE_nonredundant_vertebrates/hypergeom_per_motif.tsv \
  --cell_type ${TARGET} \
  -o /home/rshuai/research/ni-lab/kidney_finemapping/nilab/out_dir/plot_motif_enrichment
done

python3 basenji/plot/plot_allelic_imbalance_motif_enrichment.py \
  /home/rshuai/research/ni-lab/analysis/kidney_data/allelic_imbalance/motif_enrichment/tf_pseudobulk_Pseudobulk_Wilson_TF_analysis.csv \
  /home/rshuai/research/ni-lab/analysis/kidney_data/allelic_imbalance/motif_enrichment/Tubule_me_variants_neg2x_min_1e-3_diff_1_q0.01_JASPAR2020_CORE_nonredundant_vertebrates/hypergeom_per_motif.tsv \
  --cell_type Tubule \
  -o /home/rshuai/research/ni-lab/kidney_finemapping/nilab/out_dir/plot_motif_enrichment