# kidney-finemapping


## Installation
Requirements can be installed using conda.

```
conda create -n basenji_kidney_finemapping python=3.8
conda activate basenji_kidney_finemapping
conda env update --file basenji_environment.yml
pip install -e .
```
Then, to install tensorflow with GPU support, follow instructions from https://www.tensorflow.org/install/pip for your system.


## Chromatin accessibility allelic imbalance (CAAI)
### Constructing allelic imbalance and non-allelic imbalance sets
To quantify the model’s performance in assessing chromatin accessibility allelic imbalance (CAAI) at single nucleotide variants (SNVs), we first constructed a set of variants with allelic imbalance (positive set) in each of proximal tubule, loop of Henle, and distal tubule as well as a set of variants with non-allelic imbalance (negative set). To run this analysis, we used the `kidney_finemapping/basenji/make_allelic_imbalance_sets.py` script as follows:
```
TARGET=PT
NEG_MULT=7
THRESH=0.01
kidney_finemapping/basenji/make_allelic_imbalance_sets.py \
    out_dir/allelic_imbalance/data/raw/astestq10tab/all_${TARGET}q10.tsv \
    out_dir/sc_atac_seq/${TARGET}_peaks.narrowPeak \
    --neg_mult ${NEG_MULT} \
    --n_bins 20 \
    --thresh ${THRESH} \
    -o out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}
```

### Computing SAD scores for allelic imbalance and non-allelic imbalance sets
To compute SAD scores for the allelic imbalance and non-allelic imbalance sets, we used `kidney_finemapping/basenji/compute_sad.py` and merged across all chromosomes using `kidney_finemapping/basenji/merge_sad.py`.

### Allelic imbalance tasks
Two tasks were used to evaluate the model’s performance on predicting chromatin accessibility allelic imbalance (CAAI). In both tasks, we defined the model’s predictions for the reference and alternate alleles in the 192 bp bin centered at the variant as REF and ALT respectively, and we computed the predicted CAAI as REF / (REF + ALT). We tested how well predicted allelic imbalance could classify whether a variant had CAAI (i.e. discriminate CAAI variants from non-CAAI variants), as measured by both the area under the receiver operating characteristic curve (AUROC) and the area under the precision-recall curve (AUPRC). As an additional metric, among the variants in the CAAI set, we computed the AUROC for prediction of CAAI direction (REF>ALT vs ALT>REF) as measured by AUROC. The plots for these tasks can be generated using `kidney_finemapping/basenji/plot_allelic_imbalance.py`.

