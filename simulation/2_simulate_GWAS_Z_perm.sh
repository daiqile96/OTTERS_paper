#!/bin/bash
module load R/4.0.4

he2=0.1
cp=0.001
n_sim=1
n_set=1
n_thread=${NSLOTS:-1}

Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

Rfile=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_perm.R

train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

for gene in $(seq 1 ${n_row})
do

    line=`head -n $gene ${train_genes} | tail -n 1`
    ID=`echo ${line} | cut -d' ' -f 5`

    Rscript ${Rfile} ${he2} ${cp} ${n_sim} ${n_set} ${n_thread} ${ID}

done



# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N sim_Null -pe smp 1 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_perm.sh
