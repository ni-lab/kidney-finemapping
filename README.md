# kidney-finemapping

## Installation
To run the scripts in this repository, we recommend using a conda environment:
```
conda create -n basenji_kidney_finemapping python=3.8
conda activate basenji_kidney_finemapping
conda env update --file basenji_environment.yml
pip install -e .
```
Then, to install tensorflow with GPU support, follow instructions from https://www.tensorflow.org/install/pip for your system.

## Resources
TODO: Figure out how to host resources for download

All resources used in this repository can be found in the `resources` directory, and you will need to download the `resources.zip` from XXXXX. This includes the model parameters, model weights, and various data used for computing SAD scores and evaluating model performance in predicting chromatin accessibility allelic imbalance.

## SAD score analysis
### Data preprocessing
To compute SNV activity difference (SAD) scores for variants, we need the variants in VCF format. Furthermore, since we use 23 models, each holding out a different chromosome, we need to split the variants into 23 VCF files, one for each chromosome. An example of how we preprocessed the data can be found by running the `kidney_finemapping/basenji/preprocessed_finemapped_variants.py` script as follows:
```
python3 kidney_finemapping/basenji/preprocess_finemapped_variants.py \
    resources/data/220513_variants/raw/220513_gwas_replicating_loci_top_variants_susie5_finemap5_01_hg38.txt \
    --variant_set 220513 \
    -o out_dir/220513_variants/data/preprocessed
```
This script will output a VCF file for each chromosome in the `out_dir/220513_variants/data/preprocessed/snps_by_chrom` directory.

### Computing SAD scores
To compute SAD scores for variants, we used the `kidney_finemapping/basenji/compute_sad.py` script. An example of computing SAD scores for the variants on chromosome 1 is shown below:
```
CHROM=chr1
python3 kidney_finemapping/basenji/compute_sad.py \
  resources/model_params/params_sc_kidney_regression.json \
  resources/models/train_bigwigs_${CHROM}/model_best.h5 \
  out_dir/220513_variants/data/preprocessed/snps_by_chrom/${CHROM}_snps.vcf \
  -f resources/genomes/hg38.ml.fa \
  --rc \
  --shifts "1,0,-1" \
  -t resources/targets/kidney_sc_wigs_hg38.txt \
  -o out_dir/220513_variants/data/sad/${CHROM}
```
This script will output a `out_dir/220513_variants/data/sad/${CHROM}` directory containing the input VCF and an h5 file containing SAD scores for each variant. After running this script for each chromosome, we can merge the SAD scores across chromosomes using the `kidney_finemapping/basenji/merge_sad.py` script as follows:
```
python3 kidney_finemapping/basenji/merge_sad.py \
  out_dir/220513_variants/data/sad \
  --vcf \
  -o out_dir/220513_variants/data/sad/all_chrs
```

### SAD score tracks
To visualize SAD scores for a given SNV, we generated a SAD score track by plotting the model's predicted SAD score centered at every poistion for which the SNV falls within the model's receptive field. To generate these predictions, we used the `kidney_finemapping/basenji/compute_sad_shifts.py` script. An example for computing SAD score tracks for variants on chromosome 1 is shown below:
```
CHROM=chr1
python3 kidney_finemapping/basenji/compute_sad_shifts.py \
    resources/model_params/params_sc_kidney_regression.json \
    resources/models/train_bigwigs_${CHROM}/model_best.h5 \
    out_dir/220513_variants/data/preprocessed/snps_by_chrom/${CHROM}_snps.vcf \
    -f resources/genomes/hg38.ml.fa \
    --rc \
    --shifts "1,0,-1" \
    -t resources/targets/kidney_sc_wigs_hg38.txt \
    -o out_dir/220513_variants/sad_shifts/${CHROM}
```
This script will output a `out_dir/220513_variants/sad_shifts/${CHROM}` directory containing the input VCF and an h5 file containing SAD score tracks for each variant. After running this script for each chromosome, we can merge the SAD score tracks across chromosomes using the `kidney_finemapping/basenji/merge_sad_shifts.py` script as follows:
```
python3 kidney_finemapping/basenji/merge_sad.py
    out_dir/220513_variants/sad_shifts \
    -n 22 \
    --vcf \
    -o out_dir/220513_variants/sad_shifts/all_chrs
```
Then, we can plot the SAD score track for a given SNV using the `kidney_finemapping/basenji/plot_sad_shifts.py` script as follows:
```
python3 kidney_finemapping/basenji/plot/plot_sad_tracks.py \
    out_dir/220513_variants/sad_shifts/all_chrs/sad.h5 \
    -t resources/targets/kidney_sc_wigs_hg38.txt \
    --overlay \
    -o out_dir/220513_variants/plots/sad_tracks_overlay
```
This script will output a plot of the SAD score track for each SNV and save them into the `out_dir/220513_variants/plots/sad_tracks_overlay` directory.

## Chromatin accessibility allelic imbalance (CAAI)
### Constructing allelic imbalance and non-allelic imbalance sets
To quantify the model’s performance in assessing chromatin accessibility allelic imbalance (CAAI) at single nucleotide variants (SNVs), we first constructed a set of variants with allelic imbalance (positive set) in each of proximal tubule, loop of Henle, and distal tubule as well as a set of variants with non-allelic imbalance (negative set). To run this analysis, we used the `kidney_finemapping/basenji/make_allelic_imbalance_sets.py` script as follows:
```
for TARGET in PT LOH DT; do
    NEG_MULT=7
    THRESH=0.01
    python3 kidney_finemapping/basenji/make_allelic_imbalance_sets.py \
        resources/data/allelic_imbalance/raw/astestq10tab/all_${TARGET}q10.tsv \
        out_dir/sc_atac_seq/${TARGET}_peaks.narrowPeak \
        --neg_mult ${NEG_MULT} \
        --n_bins 20 \
        --thresh ${THRESH} \
        -o out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}
done
```

### Computing SAD scores for allelic imbalance and non-allelic imbalance sets
To compute SAD scores for the allelic imbalance and non-allelic imbalance sets, we used `kidney_finemapping/basenji/compute_sad.py` and merged across all chromosomes using `kidney_finemapping/basenji/merge_sad.py` as follows:

```
NEG_MULT=7
THRESH=0.01
for TARGET in PT LOH DT; do
  for SET in pos neg; do
    for CHR in {1..22}; do
      CHROM=chr${CHR}
      python3 kidney_finemapping/basenji/compute_sad.py \
        resources/model_params/params_sc_kidney_regression.json \
        resources/models/train_bigwigs_${CHROM}/model_best.h5 \
        out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}/snps_by_chrom/${SET}_set_${CHROM}.vcf \
        -f resources/genomes/hg38.ml.fa \
        --rc \
        --shifts "1,0,-1" \
        -t resources/targets/kidney_sc_wigs_hg38.txt \
        -o out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${THRESH}/${SET}_sad/${CHROM}
    done

    # Merge SAD files across chromosomes
    python3 kidney_finemapping/basenji/merge_sad.py \
      out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${THRESH}/${SET}_sad \
      --vcf \
      -o out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${THRESH}/${SET}_sad/all_chrs
  done
done
```

### Allelic imbalance tasks
Two tasks were used to evaluate the model’s performance on predicting chromatin accessibility allelic imbalance (CAAI). In both tasks, we defined the model’s predictions for the reference and alternate alleles in the 192 bp bin centered at the variant as REF and ALT respectively, and we computed the predicted CAAI as REF / (REF + ALT). We tested how well predicted allelic imbalance could classify whether a variant had CAAI (i.e. discriminate CAAI variants from non-CAAI variants), as measured by both the area under the receiver operating characteristic curve (AUROC) and the area under the precision-recall curve (AUPRC). As an additional metric, among the variants in the CAAI set, we computed the AUROC for prediction of CAAI direction (REF>ALT vs ALT>REF) as measured by AUROC. The plots for these tasks can be generated using `kidney_finemapping/basenji/plot_allelic_imbalance_tasks.py`.
```
for TARGET in PT LOH DT; do
    python3 kidney_finemapping/basenji/plot/plot_allelic_imbalance_tasks.py \
        out_dir/allelic_imbalance/sad/${TARGET}_neg7x_q0.01 \
        -t resources/targets/kidney_sc_wigs_hg38.txt \
        -o out_dir/allelic_imbalance/plots/${TARGET}_neg7x_q0.01
done
```

### Allelic imbalance motif enrichment
To examine transcription factor motifs present at sites with CAAI, we used FIMO to scan the 20 bp region surrounding each variant for matches in the JASPAR 2020 database. For each variant, we queried both the reference and alternate sequence, keeping matches with a p-value of less than 1e-3 for either the reference or alternate sequence. To identify motifs affected by variants, only matches satisfying $\left|\log_{10}\frac{p-value_{ref}}{p-value_{alt}}\ \right|\ \geq\ 1$ were analyzed.  We then computed the enrichment of motifs in the CAAI set compared to the non-CAAI set with a hypergeometric test. To do this, we used the `kidney_finemapping/basenji/allelic_imbalance_motif_enrichment.py` script as follows:
```
for TARGET in PT LOH DT; do
    python3 kidney_finemapping/basenji/allelic_imbalance_motif_enrichment.py \
        out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg7x_q0.01 \
        --motif_file resources/motif_dbs/JASPAR2020_CORE_nonredundant_vertebrates.meme \
        -o out_dir/allelic_imbalance/motif_enrichment/${TARGET}_variants_neg7x_q0.01
done
```

We then used `kidney_finemapping/basenji/plot_allelic_imbalance_motif_enrichment.py` to plot the enrichment of motifs in the CAAI set compared to the non-CAAI set. For example, to plot the enrichment of motifs in the PT CAAI set, we ran:
```
TARGET=PT
python3 kidney_finemapping/basenji/plot/plot_allelic_imbalance_motif_enrichment.py \
  resources/data/tf_pseudobulk_Pseudobulk_Wilson_TF_analysis.csv \
  out_dir/allelic_imbalance/motif_enrichment/${TARGET}_variants_neg7x_q0.01/hypergeom_per_motif.tsv \
  --cell_type ${TARGET} \
  -o out_dir/allelic_imbalance/plots/_motif_enrichment/${TARGET}_variants
```
