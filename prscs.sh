#!/bin/bash
#SBATCH -n 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-11:05                          # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                          # Memory total in MB (for all cores)
#SBATCH -o output                # File to which STDOUT will be written, including job ID
#SBATCH -e error                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0 python/2.7.12
cd /n/no_backup2/patel/yixuan/PRScs
python PRScs.py --ref_dir=/n/no_backup2/patel/yixuan/PRScs/ldblk_1kg_eur \
--bim_prefix=~/aou_bim/hapmap_aou \
--sst_file=~/sumstat.txt  \
--n_gwas=11111 \
--out_dir=~/output/ 


#score PRS using PLINK 
module load plink2

for chr in {1..22}; do
    plink2 --bfile genos/v7_chr"$chr" \
          --score output/_pst_eff_a1_b0.5_phiauto_chr"$chr".txt 2 4 6 cols=+scoresums \
          --out PRS/chr"$chr"_score
