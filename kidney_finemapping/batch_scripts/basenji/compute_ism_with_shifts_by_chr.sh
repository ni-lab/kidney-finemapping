# Run compute_ism_with_shifts.sh for each chromosome
for chr_num in 2 11 12; do
# for chr_num in {1..22}; do
    sed -i "s/CHROM=.*/CHROM=chr${chr_num}/" compute_ism_with_shifts.sh
    sbatch compute_ism_with_shifts.sh
done

# Chromosomes X and Y model
sed -i "s/CHROM=.*/CHROM=chrXY/" compute_ism_with_shifts.sh
sbatch compute_ism_with_shifts.sh
