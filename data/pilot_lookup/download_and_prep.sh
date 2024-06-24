####### Bradley June 2024 #######
# 1. Data was downloaded from google drive locally and transferred to lustre using scp

# 2. Anderson data was aggregated using the following: 
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/sc-eqtl_meta/data/pilot_lookup/Enterocytes_CD
head -n 1 cis_nominal1.cis_qtl_pairs.chr1.tsv > ../Anderson_enterocytes_sumstat.tsv
for i in {1..22}; do
    echo $i
    tail -n +2 cis_nominal1.cis_qtl_pairs.chr$i.tsv >> ../Anderson_enterocytes_sumstat.tsv
done

cp Cis_eqtls_qval.tsv ../Anderson_enterocytes_qval.tsv

# Get the sig results from Lude
cd ..
zcat franke_gut_scrna_eqtl_sumstats_enterocyte.tsv.gz | head -n 1 > franke_gut_scrna_eqtl_sumstats_enterocyte_qval.tsv
zcat franke_gut_scrna_eqtl_sumstats_enterocyte.tsv.gz | grep -v '^phenotype_id' | awk -F '\t' '$13 < 0.05' >> franke_gut_scrna_eqtl_sumstats_enterocyte_qval.tsv
gzip franke_gut_scrna_eqtl_sumstats_enterocyte_qval.tsv
