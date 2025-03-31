#######################################################
############# MTAG across traits ######################
#######################################################

formatted_dir="/path/to/formated_sumstats"
sumstat_files=$(ls $formatted_dir/*_formatted.txt | tr '\n' ',' | sed 's/,$//')

ldref_dir="/path/to/ldref_directory" #change to ancestry specific LD reference path

module load gcc/6.2.0 python/2.7.12
source py2/bin/activate
python2.7 MTAG/mtag/mtag.py \
--sumstats $sumstat_files \
--n_min 0.0 \
--out MTAG/after_stats/ \
--stream_stdout \
--ld_ref_panel $ldref_dir \
--force

#####clumping results for number of lead snps
bfilepath='/path/to/ld'
clumppath='/path/to/sumstats'
outputpath='path/of/clumped/output'

module load gcc/6.2.0 plink2/1.90
plink --bfile $bfilepath \
--clump $clumppath \
--out $outputpath \
--clump-r2 0.6 \
--clump-p1 5e-8 \
--clump-p2 0.01 \
--clump-kb 1000 \
--clump-snp-field SNP \
--clump-field mtag_pval \
