for chr_num in {1..22}
do
    sed -i "s/CHROM=.*/CHROM=chr${chr_num}/" ism_summed_pos_shifts.sh
    sbatch ism_summed_pos_shifts.sh
done

# Chromosomes X and Y model
sed -i "s/CHROM=.*/CHROM=chrXY/" ism_summed_pos_shifts.sh
sbatch ism_summed_pos_shifts.sh
