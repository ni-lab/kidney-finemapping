 for chr_num in {1..22}
do
    sed -i "s/CHROM=.*/CHROM=chr${chr_num}/" compute_sad_shifts.sh
    sbatch compute_sad.sh
done

# Chromosomes X and Y model
sed -i "s/CHROM=.*/CHROM=chrXY/" compute_sad_shifts.sh
sbatch compute_sad_shifts.sh