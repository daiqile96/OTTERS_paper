#!/bin/bash

#######################################
# SETTINGS
########################################

# # use plink to convert 
# qlogin

# load R
module load R
# load plink
module load plink/1.90b53

#########################################
# Select Genes
########################################

# Extract Gene Start and End Point from Annotation file
data_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP
anno_dir=GeneExpression/RNAseq_DLPFC/BulkDLPFC_RNAseq_2020/ROSMAP.expr.TIGAR.format_WGS_IDs.tsv
anno_file_dir=${data_dir}/${anno_dir}
sim_dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu

# Use awk to extract data
awk 'BEGIN { OFS = "\t" } { print $1, $2, $3, $4, $5 }' ${anno_file_dir} > ${sim_dir}/0_GeneLength.txt

# Use R to randomly select genes by length (select 100 genes per quantile, total 500 genes)
wk_dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu
Rscript ${wk_dir}/0_select_genes_by_length.R

#########################################
# Prepare Binary File
########################################

# train_genes=${sim_dir}/0_ReSimu_GeneAnno.txt
train_genes=${sim_dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

bin_dir=${sim_dir}/binary
mkdir ${bin_dir}
cd ${bin_dir}

test_id=${sim_dir}/test_plink_id.txt
train_id=${sim_dir}/train_plink_id.txt

for gene in $(seq 1 $n_row)
do

line=`head -n $gene ${train_genes} | tail -n 1`
ID=`echo ${line} | cut -d' ' -f 5`

# STEP 1: 
# convert the vcf genotype data for gene to bed/bim/fam file
WGS_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/VCF_per_gene
myvcf=${WGS_dir}/${ID}.vcf.gz

# make binary files from vcf file

plink --vcf ${myvcf} --make-bed --out ${ID}

# We don't have SNP ID in the vcf file so the SNPIDs are "." in the bim file.
# We use "Chrom_Pos_A1_A2" to substitute the missing SNP ID. 
cat ${ID}.bim > ${ID}_noid.bim
rm ${ID}.bim
awk 'BEGIN { OFS = "\t" } {$2=$1"_"$4"_"$5"_"$6;} 1' ${ID}_noid.bim > ${ID}.bim

# Quality control
# https://www.cog-genomics.org/plink/1.9/filter#missing
# --snps-only excludes all variants with one or more multi-character allele codes.
# --geno filters out all variants with missing call rates exceeding the provided value.
# --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01), 
# while --max-maf imposes an upper MAF bound. 
plink \
    --bfile ${ID} \
    --geno 0 --maf 0.05 --max-maf 0.45 --snps-only --make-bed \
    --out qc_${ID}

# Create VCF file from .bim, .fam and .bed file
plink --bfile qc_${ID} --recode vcf --out qc_${ID}

plink --bfile qc_${ID} --keep ${test_id} --keep-allele-order --make-bed --out test_${ID}
plink --bfile qc_${ID} --keep ${train_id} --keep-allele-order --make-bed --out train_${ID}

done 


# activate an environment to run python3
source activate demo 

# save genotype data in vcf file (qc_${ID}.vcf) to txt file (qc_${ID}_geno.txt)
for gene in $(seq 1 $n_row)
do

line=`head -n $gene ${train_genes} | tail -n 1`
ID=`echo ${line} | cut -d' ' -f 5`

python3 ${wk_dir}/0_save_geno_from_vcf.py --vcf_in=qc_${ID}.vcf --geno_out=qc_${ID}_geno.txt

done

# simulate gene expression and eQTL summary statistics
for gene in $(seq 1 $n_row)
do

line=`head -n $gene ${train_genes} | tail -n 1`
ID=`echo ${line} | cut -d' ' -f 5`

    for he2 in 0.05 0.1 ; do  
        for cp in 0.001 0.01 ; do
            Rscript ${wk_dir}/1_simulate_expression_std.R ${he2} ${cp} ${ID} 1 1
        done
    done

done



# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu -pe smp 1 -t 1 -tc 1 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/0_Simu.sh

