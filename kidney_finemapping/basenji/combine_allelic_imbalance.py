import sys
import os
import glob
import numpy as np
import scipy.stats
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

"""
Script to process allele specific counts and test allele specific imbalance

    inputs: 
        -   Directory with allele count files (vcf files generated from generate_allele_specific_counts.sh)
        -   Output directory for merged imbalance data
        -   Temporary directory storing storing individuals' heterozygous sites
    output: CSV with imbalance tested sites across heterozygous sites

    Run script like:
    python combine_allelic_imbalance.py  \
      /clusterfs/nilah/rkchung/data/atac/vcf38/ascount/ \
      /clusterfs/nilah/rkchung/data/atac/astest \
      /clusterfs/nilah/rkchung/data/atac/hetout
"""

def parse_ascounts(hetero_file):                                               
    """
    Parse allele-specific counts from a given file.

    Parameters:
    - hetero_file (str): Path to the file containing heterozygous counts.

    Returns:
    - tuple: Counts, reference allele indicator, positions, and posterior probabilities.
    """
    with open(hetero_file) as f:                                               
        # Read lines from VCF, excluding header
        het_lines = [line.strip() for line in f if not line.startswith("#")]   
        
        # Determine which lines have the reference allele first
        reffirst = np.array([1 if "0|1" in line else 0 for line in het_lines])
        
        # Extract allele-specific counts
        counts = [parse_counts(line.split("\t")[-2]) for line in het_lines]

        # Check if file has counts for two samples
        double = len(het_lines[0].split("\t")) > 11                            
        if double:
            # Get the second set of counts
            counts2 = [parse_counts(line.split("\t")[-3]) for line in het_lines]
            # Merge the two sets of counts
            counts = [[c1[0] + c2[0], c1[1] + c2[1]] for c1, c2 in zip(counts, counts2)]

        # Filter out multiallelic sites 
        valid_counts = np.array([len(c) == 2 for c in counts])                 
        counts = np.array([c for c in counts if len(c) == 2])
        reffirst = reffirst[valid_counts]

        # Extract posterior probabilities
        posterior_prob = [float(line.split("\t")[-1].split(":")[-1].split(",")[1]) for line in het_lines]
        posterior_prob = np.array(posterior_prob)[valid_counts]

        # Extract positions
        positions = [line.split("\t")[:5] for line in het_lines]
        positions = np.array(positions)[valid_counts]
        # Remove 'chr' prefix from chromosomes
        positions[:, 0] = [p.split("chr")[1] for p in positions[:, 0]]
                                                                               
    return counts, reffirst, positions, posterior_prob

def parse_counts(data):
    """
    Helper function to extract counts from a given data string.
    """
    counts_str = data.split(":")[-4].split(",")
    return [int(c) if c != "." else 0 for c in counts_str]

def compute_p_values(allele_specific_df, filter_refalt=False, count_filter=False):
    """
    Compute allelic imbalance p-values and q-values.
    Flag to remove variants with no ref and alt 

    Parameters:
    - allele_specific_df (DataFrame): DataFrame with allele-specific counts.

    Returns:
    - DataFrame: Filtered DataFrame with computed p-values.
    """
    # Convert counts to integer type
    allele_specific_df[["#Ref", "#Alt"]] = \
        allele_specific_df[["#Ref", "#Alt"]].astype(int)

    # Remove rows where the sum of reference and alternate counts is zero
    if filter_refalt:
        allele_specific_df = allele_specific_df[
            allele_specific_df["#Ref"] + allele_specific_df["#Alt"] != 0]

    # Compute fraction of reference counts to total counts
    allele_specific_df["#Ref/(#Ref+#Alt)"] = \
        allele_specific_df["#Ref"] / (allele_specific_df["#Ref"] + allele_specific_df["#Alt"])

    # Filter for sites with sufficient reads supporting heterozygous site
    if count_filter:
        allele_specific_df["Count Filter"] = \
            (allele_specific_df["#Ref"] > 1) & (allele_specific_df["#Alt"] > 1)

    # Compute binomial test p-values
    allele_specific_df["P-value"] = [
        scipy.stats.binom_test(x, n=(x+y)) if (x+y) != 0 else 1 
        for x, y in zip(allele_specific_df['#Ref'], allele_specific_df['#Alt'])]

    # Compute FDR correction
    allele_specific_df["FDR Filter"], allele_specific_df["Q-value"] = \
        fdrcorrection(allele_specific_df['P-value'], alpha=0.05)

    return allele_specific_df

def process_het_files(hetero_file, hetero_dir, out_dir, het_dir):
    """
    Process heterozygous files to extract allele-specific information and compute p-values.

    Parameters:
    - hetero_file (str): Path to the heterozygous file.
    - out_dir (str): Output directory.
    - het_dir (str): Heterozygous positions output directory.

    Returns:
    - None: Saves the processed data to a .tsv file.
    """
    print(hetero_file)
    file_name = hetero_file.split(hetero_dir)[1].split(".vcf")[0]

    # Parse counts, positions, and other data
    counts, reffirst, positions, prob = parse_ascounts(hetero_file)
    data = np.concatenate(
        (positions, prob.reshape(-1, 1), reffirst.reshape(-1, 1), counts), 
        axis=1)
    columns = ["Chr", "Pos", "-", "Ref", "Alt", "P(Het)", "Reffirst", "#Ref", "#Alt"]
    allele_specific_df = pd.DataFrame(data, columns=columns)

    # Compute p-values and other metrics
    allele_specific_df = compute_p_values(allele_specific_df, 
                                          filter_refalt=True, 
                                          count_filter=True)
    prefilt = len(allele_specific_df)

    # Heterozygous sites (for cell type specific files) called from combined cell type data
    if not "_" in "".join(file_name.split("_100bp_")):
        allele_specific_df = allele_specific_df[allele_specific_df["Count Filter"]]
        het_pos = allele_specific_df[["Chr", "Pos", "Ref", "Alt"]]
        het_pos.to_csv(f"{het_dir}/{file_name}.tsv", sep="\t", index=False)
    else:
        hetcall = "_".join(file_name.split("_")[:-1])
        if "200131" in hetcall:
            hetcall = hetcall.split("_")[0]
        het_pos = pd.read_csv(f"{het_dir}/{hetcall}q10.tsv", delim_whitespace=True)
        allele_specific_df[["Chr", "Pos"]] = allele_specific_df[["Chr", "Pos"]].astype(int)
        allele_specific_df = allele_specific_df.merge(
            het_pos, on=["Chr", "Pos", "Ref", "Alt"], how="inner")

    variant_frac = 1-(len(allele_specific_df)/prefilt)
    print(f"Finished {file_name}: {len(allele_specific_df)} variants ({variant_frac} filtered)")
    allele_specific_df.to_csv(
        os.path.join(out_dir, f"{file_name}.tsv"), sep="\t", index=False)

def merge_files(files):
    """
    Merges allelic imbalance sites across individuals. If two individuals share
    a heterozygous site, add the allele specific counts together
    
    Parameters:
    - files (list): List of files to be merged.
    
    Returns:
    - DataFrame: Merged allele-specific data.
    """
    for i, f in enumerate(files):
        df = pd.read_csv(f, sep="\t")
        
        # Add individual counts and rename columns to avoid overlap during merge
        df["IndivCount"] = 1
        df[f"P-value_{i}"] = df["P-value"]
        df[f"Q-value_{i}"] = df["Q-value"]
        df[f"#Ref/(#Ref+#Alt)_{i}"] = df["#Ref/(#Ref+#Alt)"]
        df[f"#Ref_{i}"] = df["#Ref"]
        df[f"#Alt_{i}"] = df["#Alt"]
        
        if i == 0:
            merged_df = df
        else:
            # Drop overlapping columns before merging
            cols_to_drop = ['#Alt', '#Ref', "P-value", "Q-value", "#Ref/(#Ref+#Alt)"]
            df = df.drop(cols_to_drop, axis=1)
            merged_df = merged_df.merge(df, on=["Chr", "Pos", "Ref", "Alt"], how="outer")
            
            # Combine individual counts after merging
            merged_df["IndivCount"] = merged_df["IndivCount_x"] + merged_df.get("IndivCount_y", 0)
            merged_df.drop(["IndivCount_x", "IndivCount_y"], axis=1, inplace=True)
    
    return merged_df

def process_merge(merge_indiv, cell_types, filtered, out_dir):
    """
    Process and merge individual cell type files.
    
    Parameters:
    - merge_indiv (bool): Flag to decide whether to merge or not.
    - cell_types (list): List of cell types.
    - filtered (list): List of filtered file paths.
    - out_dir (str): Output directory.
    
    Returns:
    - None: The merged files are saved to disk.
    """
    for cell_type in cell_types:
        # Filter for relevant cell type files
        if t != "combined":
            cell_files = [f for f in filtered if t in f]
        else:
            cell_files = [f for f in filtered if not any(c in f for c in cell_types[:-1])]
        
        merged_df = merge_files(cell_files)
        
        num_files = len(cell_files)
        refcols = [f"#Ref_{i}" for i in range(num_files)]
        altcols = [f"#Alt_{i}" for i in range(num_files)]
        
        # Fill NA values for counts and sum across merged files
        merged_df[refcols+altcols] = merged_df[refcols+altcols].fillna(0)
        merged_df["#Ref"] = merged_df[refcols].sum(axis=1)
        merged_df["#Alt"] = merged_df[altcols].sum(axis=1)
        merged_df.fillna({"P-value": 1, "Q-value": 1, "#Ref/(#Ref+#Alt)": 0.5}, inplace=True)
        
        # Calculate combined metrics
        merged_df = compute_p_values(merged_df,                  
                                     filter_refalt=False, 
                                     count_filter=False)

        # Save to file
        merged_df.to_csv(os.path.join(out_dir, f"all_{t}q10.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    hetero_dir = sys.argv[1]
    out_dir = sys.argv[2]
    het_dir = sys.argv[3]

    hetero_files = sorted(glob.glob(os.path.join(hetero_dir, "*q10.vcf")), key=lambda x: len(x.split("_")))

    print(hetero_files)
    generate_het = True
    if generate_het: # generates allele counts for heterozygous sites (from vcf)
        for hetero_file in hetero_files:
            process_het_files(hetero_file, hetero_dir, out_dir, het_dir)

    merge_indiv = True
    filtered = glob.glob(out_dir+"*q10.tsv")
    cell_types = list(filter(lambda c: len(c)<7, set(f.split("q10.tsv")[0].split("_")[-1] for f in filtered)))
    filtered = list(filter(lambda c: len(c.split("_100bp_"))!=2, filtered))
    filtered = list(filter(lambda c: len(c.split("all_"))!=2, filtered))
    print(cell_types)
    print(filtered)
    cell_types += ["combined"]
    if merge_indiv:
        process_merge(cell_types, filtered, out_dir)
