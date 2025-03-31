#!/bin/bash
#SBATCH -n 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-11:05                          # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                          # Memory total in MB (for all cores)
#SBATCH -o output                # File to which STDOUT will be written, including job ID
#SBATCH -e error                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

module load gcc/6.2.0 python/3.6.0
cd /n/no_backup2/patel/yixuan/PRScsx
python3 PRScsx.py --ref_dir=~/PRScs/ \
--bim_prefix=~/aou_bim/hapmap_aou \
--sst_file=~/pheno_AFR.txt,\
~/pheno_AMR.txt,\
~/pheno_EAS.txt,\
~/pheno_EUR.txt \
--n_gwas=11111,22222,33333,44444 \
--pop=AFR,AMR,EUR,EAS \
--out_dir=~/output/ \
--out_name=pheno \
--chrom=${1} \
--meta=True

#score PRS using PLINK 
module load plink2

for anc in AFR AMR EAS EUR META; do               
  for chr in {1..22}; do
    plink2 --bfile genos/v7_chr"$chr" \
          --score output/pheno_"$anc"_pst_eff_a1_b0.5_phiauto_chr"$chr".txt 2 4 6 cols=+scoresums \
          --out PRS/"$anc"/chr"$chr"_score
  done
done
