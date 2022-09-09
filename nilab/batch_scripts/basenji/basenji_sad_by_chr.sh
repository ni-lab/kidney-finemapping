 for chr_num in {1..22}
do
    sed -i "s/CHROM=.*/CHROM=chr${chr_num}/" basenji_sad.sh
    sbatch basenji_sad.sh
done

# Chromosomes X and Y model
sed -i "s/CHROM=.*/CHROM=chrXY/" basenji_sad.sh
sbatch basenji_sad.sh