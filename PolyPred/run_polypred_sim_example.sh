#!/bin/bash

# load modules and conda environment
module load anaconda3/2022.05
conda init bash
conda activate polyfun

polyfun_path="path_to_polyfun_function"
plink="path_to_plink" 

# Specify simulated data
sim_path=$(cat sim_path_${sim_index}.txt)

rep_index=${SGE_TASK_ID}

for trait_index in 4;
do
#created munged sumstats
python $polyfun_path/munge_polyfun_sumstats.py \
    --sumstats $sim_path/rep${rep_index}_${trait_index}.txt \
    --out $sim_path/rep${rep_index}_${trait_index}.sumstats.munged.parquet \
    --n 337491

python $polyfun_path/create_finemapper_jobs.py \
    --sumstats $sim_path/rep${rep_index}_${trait_index}.sumstats.munged.parquet \
    --n 337491 \
    --method susie \
    --non-funct \
    --max-num-causal 5 \
    --allow-missing \
    --out-prefix $sim_path/rep${rep_index}_${trait_index}_polyfun_output \
    --jobs-file $sim_path/rep${rep_index}_${trait_index}_jobs.txt

cd $sim_path
R --vanilla --args $sim_path ${rep_index} ${trait_index}<< "REND"
  args = commandArgs(trailingOnly = TRUE)

library(data.table)
jobs=fread(paste0(args[1],"/rep",args[2],"_",args[3],"_jobs.txt"),header=F,sep="\t")
path_ld="path_to_downloaded_UKBB_LD"
jobs.new=gsub("https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/",path_ld,jobs$V1)
fwrite(as.data.frame(jobs.new),file=paste0(args[1],"/rep",args[2],"_",args[3],"_jobs_new.txt"),row.names=F, col.names=F)
REND


#run each and every fine-mapping job (could take several hours or more)
module load R/3.6.3
bash $sim_path/rep${rep_index}_${trait_index}_jobs_new.txt
python $polyfun_path/aggregate_finemapper_results.py \
    --out-prefix $sim_path/rep${rep_index}_${trait_index}_polyfun_output \
    --sumstats $sim_path/rep${rep_index}_${trait_index}.sumstats.munged.parquet \
    --out $sim_path/rep${rep_index}_${trait_index}_polyfun_output.agg.txt.gz \
    --adjust-beta-freq \
    --allow-missing-jobs

python $polyfun_path/polypred.py --combine-betas --betas $sim_path/rep${rep_index}_${trait_index}.txt,$sim_path/rep${rep_index}_${trait_index}_polyfun_output.agg.txt.gz --pheno /SFS/project/genetics/guobin/scratch/mtPRS_finemap/pgx_info/target.train.pheno.txt --output-prefix $sim_path/rep${rep_index}_${trait_index}_polypred --plink-exe $plink /SFS/project/genetics/guobin/scratch/mtPRS_finemap/pgx_info/target.train.plink.bed

python $polyfun_path/polypred.py \
    --predict \
    --betas $sim_path/rep${rep_index}_${trait_index}_polypred.betas \
    --output-prefix $sim_path/rep${rep_index}_${trait_index}_polypred.predictions \
    --plink-exe $plink \
    target.plink.bed

done

rm *log
rm *gz












