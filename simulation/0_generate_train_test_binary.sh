#!/bin/bash

# Load Plink
module load plink/1.90b53

# Set Data directory
sim_dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu

# Generate sample IDs
test_id=${sim_dir}/test_plink_id.txt
train_id=${sim_dir}/train_plink_id.txt
awk 'BEGIN { OFS = "\t" } { print $1, $1 }' ${sim_dir}/test_id.txt > ${test_id}
awk 'BEGIN { OFS = "\t" } { print $1, $1 }' ${sim_dir}/train_id.txt > ${train_id}

# train_genes=${sim_dir}/0_ReSimu_GeneAnno.txt
train_genes=${sim_dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

bin_dir=${sim_dir}/binary

cd ${bin_dir}

for gene in $(seq 1 $n_row)
do

line=`head -n $gene ${train_genes} | tail -n 1`
ID=`echo ${line} | cut -d' ' -f 5`


echo ${ID}
# generate test and train samples
plink --bfile qc_${ID} --keep ${test_id} --keep-allele-order --make-bed --out test_${ID}
plink --bfile qc_${ID} --keep ${train_id} --keep-allele-order --make-bed --out train_${ID}

done 

