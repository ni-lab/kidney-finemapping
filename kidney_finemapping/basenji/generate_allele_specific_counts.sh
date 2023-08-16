#!/bin/bash
# Gets allele specific count from bam file for heterozygous sites (in inputed genotype) 
# Takes four arguments:
#   1. Data directory containing reference genome fasta and imputed genotype vcfs 
#   2. Output directory outputs (results saved to outdir/het)
#   3. Bam file directory containing ATAC seq bams  
#   4. Names of each individual that has been genotyped
#   5. Name of bam file directories that matches the order of the individuals in (4)
# Example:
#   bash generate_allele_specific_counts.sh /clusterfs/nilah/rkchung/data/atac/vcf38/ \
#         /clusterfs/nilah/rkchung/data/atac/vcf38/ascount/ \
#         /clusterfs/nilah/kidney_ATACseq/WASP/ \
#         "204686210102_R01C01 204686210102_R02C01 204686210102_R03C01 204686210102_R03C01" \
#         "200131cortex 200317combined 200707_100bp_combined 200707combined"
datadir=$1
outdir=$2
bamdir=$3
indiv=($4)
experdir=($5)

mkdir $outdir"mp"
mkdir $outdir"het"

echo "indivs ${indiv[@]}"
echo "expers ${experdir[@]}"

# Forms allele specific counts for heterozygous sites using ATAC seq bam files and imputed genotypes
for i in "${!indiv[@]}"; do
    echo "Current individual ${indiv[$i]}"
    dir=$datadir"${indiv[$i]}"/

    # Get all bam files fron experiment directory
    bs=($(ls $bamdir${experdir[$i]}/*.bam))
    #bams=($(printf '%s\n' "${bs[@]}" | tac | tr '\n' ' ')) # reverse list
    bams=(${bs[@]})
    for j in "${!bams[@]}"; do # iterate through all bams
        # Gets short name experiment name + cell type 
        short=$(echo ${bams[$j]} | sed -n "s/.*sorted\([^']*\)\.bam.*/\1/p")
        short=${experdir[$i]}$short
        echo "Experiment name: "$short

        short=$short"q10_test"
        # Generate vcf if not created
        if [ ! -f $outdir"het/"$short".vcf" ]; then
            echo "Generating vcf using bcftools mpileup and call"
            bcftools mpileup -Ou -d 1000 -q 10 -a FORMAT/AD,FORMAT/DP -f $datadir"hg38.fa" ${bams[$j]} | bcftools call -m > $outdir"mp/"$short".vcf"
            bgzip $outdir"mp/"$short".vcf"
            bcftools index -f -c $outdir"mp/"$short".vcf.gz" 

            echo "Filtering out homozygous sites from imputed genotype"
            bcftools view -h $dir"merged_lifted.vcf.gz" > $dir"merged_lifted_het.vcf"
            bcftools view $dir"merged_lifted.vcf.gz" | grep "1|0\|0|1" >> $dir"merged_lifted_het.vcf"
            bgzip -f $dir"merged_lifted_het.vcf"
            bcftools index -f -c $dir"merged_lifted_het.vcf.gz" 

            echo "Merge atac-seq mpileup and heterozygous genotyped sites"
            bcftools isec -c all -p $outdir"mp/"$short -Oz $outdir"mp/"$short".vcf.gz" $dir"merged_lifted_het.vcf.gz"
            bcftools merge --merge all $outdir"mp/"$short/0002.vcf.gz $outdir"mp/"$short/0003.vcf.gz > $outdir"het/"$short".vcf" 

            echo "Cleaning up bcftools isec"
            rm -rf $outdir"mp/"$short
        fi
    done
done
