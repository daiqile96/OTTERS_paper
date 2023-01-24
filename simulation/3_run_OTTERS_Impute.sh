#!/bin/bash

# activate the environment
conda activate otters

# Load R to perform lassosum
module load R/4.0.4

# Load Plink
module load plink/1.90b53

# Load Tabix
module load tabix

N_THREADS=1
OTTERS_DIR=/home/qdai/YangFSSdata2/qdai/tool/OTTERS
Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

for he2 in 0.01 0.05 0.1 ; do  
    
    for cp in 0.001 0.01 ; do

		# Start to run OTTERS
		Result_Dir=${Simu_Dir}/Results/std_cp_${cp}_He2_${he2}

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
			# ID to calculate GReX
			sample_dir=${Simu_Dir}/test_plink_id.txt
			# Output dir
			grex_dir=${Result_Dir}/GReX/${ID}
			mkdir -p ${grex_dir}

		    python3 ${OTTERS_DIR}/imputing.py \
				--OTTERS_dir=${OTTERS_DIR} \
				--weight_dir=${Result_Dir}/${ID} \
				--models=SDPR,PRScs,lassosum,P0.001,P0.05,FUSION \
				--anno_dir=${exp_anno} \
				--geno_dir=${geno_dir} \
				--out_dir=${grex_dir} \
				--chrom=${CHROM} \
				--samples=${sample_dir} \
				--thread=${N_THREADS}

		done

	done

done

# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu_StageII -pe smp 3 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/3_run_OTTERS_Impute.sh


# for he2 in 0.05 0.1 ; do  
    
#     for cp in 0.001 0.01 ; do

#     	Rscript /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/4_evaluate_R2.R ${he2} ${cp}

#     done

#  done
