# kidney-finemapping


## Installation
Requirements can be installed using conda.

```
conda env create -f environment.yml
conda install tensorflow-gpu
pip install -e .
```

## Allelic imbalance task
To quantify the model’s performance in assessing chromatin accessibility allelic imbalance (CAAI) at single nucleotide variants (SNVs), we constructed a set of variants with allelic imbalance in each of proximal tubule, loop of Henle, and distal tubule. To run this analysis, `kidney_finemapping/basenji/make_allelic_imbalance_sets.py`