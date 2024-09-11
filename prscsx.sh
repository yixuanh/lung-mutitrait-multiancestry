SUMSTAT_dir="/path/to/sumstats"
sumstat_files=$(ls $formatted_dir/*_formatted.txt | tr '\n' ',' | sed 's/,$//')

ldref_dir="/path/to/ldref_directory" #change to ancestry specific LD reference path

#get weights using PRScsx
python3 PRScsx.py \
               --ref_dir=/ldrefs \
               --bim_prefix=/hapmap_aou \
               --sst_file=$SUMSTAT_dir \
               --n_gwas=32658,18173,341204,1376071 \
               --pop=AFR,AMR,EAS,EUR \
               --out_dir=/mnt/data/output/ \
               --meta=True \
               --out_name=  
               
#score PRS using PLINK 
for anc in AFR AMR EAS EUR; do               
  for chr in {1..22}; do
    plink2 --bfile genos/aou_v7_chr"$chr" \
          --score SNPweights/PRSCSX/"$chr"/"$anc"/_pst_eff_a1_b0.5_phiauto_chr"$chr".txt 2 4 6 cols=+scoresums \
          --out PRS/V3/Asthma/PRSCSX/"$anc"/chr"$chr"_score_"$anc"
  done
done
