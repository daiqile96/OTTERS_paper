#!/bin/bash

# activate the environment
conda activate otters

# Load R to perform lassosum
module load R/4.0.4

# Load Plink
module load plink/1.90b53

# Load Tabix
module load tabix

# GET JOB INFO FROM LIST OF TRAIN JOBS 
i=${SGE_TASK_ID}

# he2=0.01
# for cp in 0.001 0.01 ; do
#     for N in 50 75 100 150 200; do 
#         echo -e "${he2}\t${cp}\t${N}" >> ${train_jobs_path}
#     done 
# done

# he2=0.01
# for cp in 0.001 0.01 ; do
#     for N in 300 400 500; do 
#         echo -e "${he2}\t${cp}\t${N}" >> ${train_jobs_path}
#     done 
# done

# he2=0.05
# for cp in 0.001 0.01 ; do
#     for N in 25 50 75 100 ; do 
#         echo -e "${he2}\t${cp}\t${N}" >> ${train_jobs_path}
#     done 
# done

# he2=0.1
# for cp in 0.001 0.01 ; do
#     for N in 10 25 50 75 ; do 
#         echo -e "${he2}\t${cp}\t${N}" >> ${train_jobs_path}
#     done 
# done

#train_jobs_path=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/3_new_simulation/4_sim_GWAS_Zscore_job.txt
train_jobs_path=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_allN.txt

line=`head -n $i ${train_jobs_path} | tail -n1`
he2=`echo $line | cut -d' ' -f 1`
cp=`echo $line | cut -d' ' -f 2`
N=`echo $line | cut -d' ' -f 3`
N_THREADS=1

# set up my OTTERS directory directory
OTTERS_DIR=/home/qdai/YangFSSdata2/qdai/tool/OTTERS
Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/


# Start to run OTTERS
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

Result_Dir=${Simu_Dir}/Results/std_cp_${cp}_He2_${he2}
Sub_Result_Dir=${Result_Dir}/TWAS/N${N}

for gene in $(seq 1 ${n_row})
    
do

	line=`head -n $gene ${train_genes} | tail -n 1`
	ID=`echo ${line} | cut -d' ' -f 5`
	CHROM=`echo ${line} | cut -d' ' -f 1`

	echo ${ID}

	# Input for OTTERS STAGE II Imputation 
	# Annotation File 
	exp_anno=${Result_Dir}/${ID}/${ID}_anno.txt
	# Genotype data from LD reference panel
	geno_dir=${Simu_Dir}/binary/test_${ID}

	for set in $(seq 1 10)
	do 
			# Set output directory for STAGE II
			sub_twas_dir=${Sub_Result_Dir}/${ID}_Set${set}
			mkdir -p ${sub_twas_dir}

			# Input for OTTERS STAGE II
			# GWAS summary statistics 
			gwas_sst_file=${Simu_Dir}/Training/std_cp_${cp}_He2_${he2}/${ID}_Hp0.025_Set${set}_N${N}K_GWAS_Zscore_1.txt
			

			python3 ${OTTERS_DIR}/testing.py \
				--OTTERS_dir=${OTTERS_DIR} \
				--weight_dir=${Result_Dir}/${ID} \
				--models=FUSION \
				--anno_dir=${exp_anno} \
				--geno_dir=${geno_dir} \
				--out_dir=${sub_twas_dir} \
				--gwas_file=${gwas_sst_file} \
				--chrom=${CHROM} \
				--thread=$N_THREADS
	done

done

# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu_StageII -pe smp 3 -t 1-32 -tc 16 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/3_run_OTTERS_StageII_FUSION.sh
# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu_StageII -pe smp 1 -t 1-10 -tc 10 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/3_run_OTTERS_StageII.sh


# for i in $(seq 1 16)
# do 

# 	train_jobs_path=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_allN.txt

# 	line=`head -n $i ${train_jobs_path} | tail -n1`
# 	he2=`echo $line | cut -d' ' -f 1`
# 	cp=`echo $line | cut -d' ' -f 2`
# 	N=`echo $line | cut -d' ' -f 3`

# 	Rscript /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/4_evaluate_Power.R ${he2} ${cp} ${N}

# done
