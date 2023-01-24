# activate the environment
conda activate otters

# Load R to perform lassosum
module load R/4.0.4

# Load Plink
module load plink/1.90b53

# Load Tabix
module load tabix

Batch=$SGE_TASK_ID
# set number of threads to be used
N_THREADS=1

# set up my OTTERS directory and SDPR directory
OTTERS_DIR=/home/qdai/YangFSSdata2/qdai/tool/OTTERS
SDPR_DIR=/home/qdai/YangFSSdata2/qdai/tool/SDPR2

# make sure the dynamic libraries of SDPR are not changed (For SDPR)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib

# prevent automatically using  all available cores on a compute node (For SDPR and PRS-CS)
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS


# for he2 in 0.05 0.1 ; do  

for he2 in 0.01 ; do  
	
    for cp in 0.001 0.01 ; do

	# Start to run OTTERS
	Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/
	Result_Dir=${Simu_Dir}/Results/std_cp_${cp}_He2_${he2}
	mkdir -p ${Result_Dir}

	# train_genes=${Simu_Dir}/0_ReSimu_GeneAnno.txt
	# we increase the number of genes here
	train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt

	n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

		for gene in $(seq $(( 10*Batch - 9 )) $(( 10*Batch )))
		do


		line=`head -n $gene ${train_genes} | tail -n 1`
		ID=`echo ${line} | cut -d' ' -f 5`
		CHROM=`echo ${line} | cut -d' ' -f 1`
		START=`echo ${line} | cut -d' ' -f 2`
		END=`echo ${line} | cut -d' ' -f 3`

		echo ${ID}

		Sub_Result_Dir=${Result_Dir}/${ID}
		mkdir -p ${Sub_Result_Dir}
		echo -e "CHROM\tGeneStart\tGeneEnd\tTargetID" > ${Sub_Result_Dir}/${ID}_anno.txt
		echo -e "${CHROM}\t${START}\t${END}\t${ID}" >> ${Sub_Result_Dir}/${ID}_anno.txt

		# Input for OTTERS STAGE I 
		# Annotation File 
		exp_anno=${Sub_Result_Dir}/${ID}_anno.txt
		# Genotype data from LD reference panel
		geno_dir=${Simu_Dir}/binary/train_${ID}
		# eQTL summary statistics 
		sst_file=${Simu_Dir}/Training/std_cp_${cp}_He2_${he2}/${ID}_sst_1.txt
		# Input for OTTERS STAGE II
		# GWAS summary statistics 
		# gwas_sst_file=Exp_GWASSumStats.txt

		# Set LD-clumping threshold in STAGE I
		clump_r2=0.99
		# Set output directory for STAGE I
		out_dir=${Sub_Result_Dir}

		# STAGE I
		# train eQTL weights using P+T, lassosum, SDPR and PRS-CS. 
		# It may take several minutes to complete.
		python3 ${OTTERS_DIR}/training.py \
		--OTTERS_dir=${OTTERS_DIR} \
		--SDPR_dir=${SDPR_DIR} \
		--anno_dir=${exp_anno} \
		--geno_dir=${geno_dir} \
		--sst_file=${sst_file} \
		--out_dir=${out_dir} \
		--chrom=${CHROM} \
		--r2=${clump_r2} \
		--models=PT,lassosum,SDPR,PRScs \
		--thread=$N_THREADS

		done

	done

done


# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu_OTTERS -pe smp 2 -t 1-50 -tc 50 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/1_run_OTTERS.sh
