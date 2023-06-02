# kidney-finemapping


## Installation
Requirements can be installed using conda.

```
conda env create -f environment.yml
conda install tensorflow-gpu
pip install -e .
```

## Allelic imbalance task

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

