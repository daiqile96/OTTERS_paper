
# transfer summary statistics
rsync -avzh qdai@hgcc.genetics.emory.edu:/mnt/YangFSS/data2/GTEx/GTEx_V8/ExpressionFiles/eQTL/GTEx_Analysis_v8_eQTL_all_associations/Whole_Blood.allpairs.txt.gz /home/qdai8/projects/OTTERS/Data/ReRDA/

# Format raw eQTL data
conda activate otters
data_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/
work_dir=/home/qdai8/projects/OTTERS/Scripts/3_ReRDA
out_unsort=$data_dir/eQTL/GTEx_V8_Whole_Blood/Whole_Blood.allpairs.txt

python3 $work_dir/0.1_covert_z.py \
--stat_dir=$data_dir/Whole_Blood.allpairs.txt.gz \
--out_dir=${out_unsort} \
--N=670

# sort eQTL data by CHROM and POS
nthread=4
buffersize=8G
out_file=$data_dir/eQTL/GTEx_V8_Whole_Blood/Whole_Blood.allpairs.sort.txt
sort -nk1 -nk2 --parallel ${nthread} -S ${buffersize} ${out_unsort} >> ${out_file} && rm ${out_unsort} || error_exit 'Unable to sort filtered output file.'