
# create file for each chromosome 
for chr in {1..22}
do
	for method in lassosum P0.001 P0.05 PRScs SDPR FUSION
	do 
		out_twas=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/${method}.txt
		echo -e 'CHROM\tGeneStart\tGeneEnd\tTargetID\tn_snps\tFUSION_Z\tFUSION_PVAL' > ${out_twas}
	done
done

# merge results from different batches
train=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/train_group.txt
for set in {1..207}
do
	line=`head -n $set ${train} | tail -n 1`
	chr=`echo ${line} | cut -d' ' -f 1`
	batch=`echo ${line} | cut -d' ' -f 2`
	
	twas_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/Batch${batch}/TWAS

	for method in lassosum P0.001 P0.05 PRScs SDPR FUSION
	do 
		out_twas=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/${method}.txt
		tmp_twas=${twas_dir}/${method}.txt
		tail -n +2 ${tmp_twas} >> ${out_twas}
	done
done


# # create file for each chromosome 
# for chr in {1..22}
# do
# bmi_dir=/home/qdai8/projects/OTTERS/Data/RDA/Results/CHR${chr}/BMI
# mkdir -p ${bmi_dir}
# 	for method in lassosum P0.001 P0.05 PRScs SDPR 
# 	do 
# 		out_twas=${bmi_dir}/${method}.txt
# 		echo -e 'CHROM\tGeneStart\tGeneEnd\tTargetID\tn_snps\tFUSION_Z\tFUSION_PVAL' > ${out_twas}
# 	done
# done

# # merge results from different batches
# train=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/train_group.txt
# for set in {1..207}
# do
# 	line=`head -n $set ${train} | tail -n 1`
# 	chr=`echo ${line} | cut -d' ' -f 1`
# 	batch=`echo ${line} | cut -d' ' -f 2`
# 	bmi_dir=/home/qdai8/projects/OTTERS/Data/RDA/Results/CHR${chr}/BMI
# 	twas_dir=/home/qdai8/projects/OTTERS/Data/RDA/Results/CHR${chr}/Batch${batch}/TWAS_BMI

# 	for method in lassosum P0.001 P0.05 PRScs SDPR 
# 	do 
# 		out_twas=${bmi_dir}/${method}.txt
# 		tmp_twas=${twas_dir}/${method}.txt
# 		tail -n +2 ${tmp_twas} >> ${out_twas}
# 	done
# done