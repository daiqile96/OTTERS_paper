#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=40:00:00
#SBATCH --partition=week-long-cpu,yanglab,day-long-cpu
#SBATCH --job-name=FUSION
#SBATCH --error=/home/qdai8/projects/logs/job.%J.error
#SBATCH --output=/home/qdai8/projects/logs/job.%J.out

# Load modules 
module purge
module load R/4.0.3

Simu_Dir=/home/qdai8/projects/OTTERS/Data/ReSimu/
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt

for gene in {1..500} ; do
	
	line=`head -n $gene ${train_genes} | tail -n 1`
	ID=`echo ${line} | cut -d' ' -f 5`

	for cp in 0.001 0.01 ; do

	    for he2 in 0.01 0.05 0.1 ; do 

	    Rscript /home/qdai8/projects/OTTERS/Scripts/2_ReSimu/0_pre_fusion.R ${cp} ${he2} ${ID}

	    cd /home/qdai8/projects/OTTERS/Data/ReSimu/binary

	    b_dir=train_${ID}
		tmp_dir=${ID}_${cp}_${he2}.tmp
		o_dir=${ID}_${cp}_${he2}

		Rscript /home/qdai8/projects/OTTERS/Tool/fusion_twas-master/FUSION.compute_weights.R \
		--bfile $b_dir \
		--tmp $tmp_dir \
		--out $o_dir \
		--hsq_p 1 \
		--PATH_plink /home/qdai8/projects/OTTERS/Tool/plink \
		--PATH_gcta /home/qdai8/projects/OTTERS/Tool/fusion_twas-master/gcta_nr_robust \
		--PATH_gemma /home/qdai8/projects/OTTERS/Tool/gemma-linux \
		--models top1,lasso,enet,blup

	    done 

	done

done


