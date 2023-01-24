# Load R 
module load R/4.0.4


Simu_Dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/
train_genes=${Simu_Dir}/0_ReSimu_GeneAnno_100_per_group.txt
n_row=`wc -l ${train_genes} | cut -d' ' -f 1`

Batch=$SGE_TASK_ID

for gene in $(seq $(( 100*Batch - 99 )) $(( 100*Batch )))

	do

	line=`head -n $gene ${train_genes} | tail -n 1`
	ID=`echo ${line} | cut -d' ' -f 5`

	for he2 in 0.01 0.05 0.1 ; do  

        for cp in 0.001 0.01 ; do

            Rscript /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/6_save_fusion.R ${cp} ${he2} ${ID}

        done

     done

done

# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N ReSimu_OTTERS -pe smp 2 -t 1-5 -tc 5 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/3_save_fusion.sh
