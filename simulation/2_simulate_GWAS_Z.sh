#!/bin/bash
hp2=0.025
module load R/4.0.4

# GET JOB INFO FROM LIST OF TRAIN JOBS 
i=${SGE_TASK_ID}

train_jobs_path=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_allN_with_Batch.txt

# he2=0.01
# for cp in 0.001 0.01 ; do
#     for N in 50 75 100 150 200 300 400 500; do 
#         for Batch in {1..50}; do
#            echo -e "${he2}\t${cp}\t${N}\t${Batch}" >> ${train_jobs_path}
#         done
#     done 
# done

# he2=0.05
# for cp in 0.001 0.01 ; do
#     for N in 25 50 75 100; do 
#         for Batch in {1..50}; do
#            echo -e "${he2}\t${cp}\t${N}\t${Batch}" >> ${train_jobs_path}
#         done
#     done 
# done

# he2=0.1
# for cp in 0.001 0.01 ; do
#     for N in 10 20 30 40 ; do 
#         for Batch in {1..50}; do
#            echo -e "${he2}\t${cp}\t${N}\t${Batch}" >> ${train_jobs_path}
#         done
#     done 
# done


line=`head -n $i ${train_jobs_path} | tail -n1`
he2=`echo $line | cut -d' ' -f 1`
cp=`echo $line | cut -d' ' -f 2`
N=`echo $line | cut -d' ' -f 3`
Batch=`echo $line | cut -d' ' -f 4`
n_sim=1
n_set=10
n_thread=${NSLOTS:-1}


Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

Rfile=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z.R

for gene in $(seq $(( 10*Batch - 9 )) $(( 10*Batch )))
    
    do

    line=`head -n $gene ${train_genes} | tail -n 1`
    ID=`echo ${line} | cut -d' ' -f 5`

    Rscript ${Rfile} ${he2} ${cp} ${N} ${hp2} ${n_sim} ${n_set} ${n_thread} ${ID}

done



# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N sim_Zscore -pe smp 2 -t 1-1600 -tc 16 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z.sh
