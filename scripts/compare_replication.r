########################################################################################
########### Bradley - June 2024. Checking replication of eQTLs across cohorts #########
########################################################################################

# Set up
library(ggplot2)
library(dplyr)
library(ggvenn)
library(ComplexUpset)
library(tidyr)
library(qvalue)
library(remotes)
library(purrr)
library(ComplexHeatmap)
setwd("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/sc-eqtl_meta/")
pathOut = "results/pilot_lookup"

# Define file names and directory
fdir = "data/pilot_lookup"
load_info=list()
load_info["Georges_enterocyte_qval"] = "20240513_Liege_scrna_eqtl_gene_sumstats_enterocyte.tsv"
load_info["Georges_enterocyte_tier_2_qval"] = "20240516_Liege_scrna_eqtl_module_sumstats_enterocyte_tier_2.tsv"
load_info["Anderson_enterocyte_qval"] = "Anderson_enterocytes_qval.tsv"
#load_info["Anderson_enterocyte_sumstat"] = "Anderson_enterocytes_sumstat.tsv"
load_info["Franke_enterocyte_qval"] = "franke_gut_scrna_eqtl_sigpairs_enterocyte.tsv.gz"
load_info["Franke_enterocyte_sumstat_qval"] = "franke_gut_scrna_eqtl_sumstats_enterocyte.tsv.gz" 

# Load in
all = vector("list", length = length(load_info))
for (i in 1:length(all)){
    f = load_info[[i]]
    fname = names(load_info)[i]
    print(fname)
    if (grepl(".gz", f)){ # Unzip if required
        if (grepl("franke",f) & grepl("sigpairs", f)){
            all[[i]] = read.csv(paste0(fdir, "/", f), sep = "\t", header=F)
            colnames(all[[i]]) = c("variant_id", "phenotype_id") # rename cols if Franke sig pairs
        } else {
            all[[i]] = read.csv(gzfile(paste0(fdir, "/", f)), sep = "\t")
        }
    } else {
        all[[i]] = read.csv(paste0(fdir, "/", f), sep = "\t")
    }
    all[[i]]$dataset = fname # Add dataset name
    names(all)[i] = fname
    # Rename columns
    if("gene_id" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "gene_id")] = "phenotype_id"
    }
    if("phe_id" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "phe_id")] = "phenotype_id"
    }
    if("qval" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "qval")] = "qvalue"
    }
    if("pval_beta" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "pval_beta")] = "p_beta"
    }
    if("adj_beta_pval" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "adj_beta_pval")] = "p_beta"
    }
    if("pval_perm" %in% colnames(all[[i]])){
        colnames(all[[i]])[which(colnames(all[[i]]) == "pval_perm")] = "p_beta" # CHECK THIS (APPLIES TO FRANKE DATASET OF SIG QVALS)
    }
}

################## 1. Look at the overlap of significant eGenes, irrespective of direction ###################
anderson_sig = all[["Anderson_enterocyte_qval"]][all[["Anderson_enterocyte_qval"]]$qval < 0.05,]
lude_sig = 
list_data <- list(Georges = unique(all$Georges_enterocyte_qval$gene_id), 
                    Georges_tier2 = unique(all$Georges_enterocyte_tier_2_qval$phe_id), 
                    Anderson = unique(anderson_sig$phenotype_id), Franke = unique(all$Franke_enterocyte_qval$phenotype_id))

p = ggvenn(list_data)
ggsave(filename = paste0(pathOut, "/venn_qval_hit_replication.png"), plot = p, width = 8, height = 6)

# Also make an upset plot of this
unique_ids <- unique(unlist(list_data))
df <- sapply(list_data, function(set) unique_ids %in% set)
df <- as.data.frame(df)
rownames(df) <- unique_ids
df_long <- df %>%
  mutate(ID = rownames(df)) %>%
  pivot_longer(-ID, names_to = "Set", values_to = "Present")

# Create the upset plot using ComplexUpset
up = upset(df, colnames(df), encode_sets=T)
ggsave(filename = paste0(pathOut, "/upset_qval_hit_replication.png"), plot = up, width = 8, height = 6)


upset_plot <- upset(
  df_long, 
  intersect = list_data, 
  name = 'Set', 
  width_ratio = 0.1, 
  annotations = list(
    'ID' = ggplot() + geom_bar(aes(x = Present))
  )
)

################## 3. Calculate the replication of hit from one dataset into another ###################
# ref_page = https://rdrr.io/github/kauralasoo/seqUtils/src/R/qtl_postProcess.R
#' Currently expects the following columns: gene_id, qvalue and p_beta.
#' @param table1 First table maximum p-values per feature.
#' @param table2 Second table of maximum p-values per feature.
#' @param qvalue_thresh qvalue threshold for table1.
#' @return None
#' @author Kaur Alasoo
#' @export
calculatePi1 <- function(table1, table2, qvalue_thresh = 0.1, feature_id = "gene_id"){
  #Identify significant hits from first table
  table1_hits = dplyr::filter(table1, qvalue < qvalue_thresh)
  #Extract the same genes from the second table
  table2_hits = dplyr::semi_join(table2, table1_hits, by = feature_id)
  #Estimate the proportion of replicated qtls
  pi1 = 1 - qvalue::qvalue(table2_hits$p_beta)$pi0
  return(pi1)
}

# Subset input for those with required columns
req_cols = c("phenotype_id", "qvalue", "p_beta")

contains_all_req_cols <- function(df, req_cols) {
  all(req_cols %in% names(df))
}

all_use <- keep(list_of_dfs, ~ contains_all_req_cols(.x, req_cols))


# generate a pairwise matrix for replication in one dataset vs another
pigrid = matrix(length(all_use), length(all_use))
rownames(pigrid) = names(all_use)
colnames(pigrid) = names(all_use)

for ( r in 1:length(all_use) ){
    for ( c in 1:length(all_use) ){
        if(c != r){
            # Ref dataset is the row
            print(paste0("Ref ", names(all_use)[r], " vs query ", names(all_use)[c]))
            pi1 = calculatePi1(table1 = all[[r]][,req_cols],
                                    table2 = all[[c]][,req_cols],
                                    qvalue_thresh = 0.05,
                                    feature_id = "phenotype_id")
            # Add
            pigrid[r,c] = pi1
        }
    }
}

# Plot this:
library(ComplexHeatmap)
Heatmap(mat_df,
        name = "Values",  # Name of the matrix
        col = colorRamp2(c(0, 1), c("white", "blue")),  # Color scheme (adjust as needed)
        show_row_names = TRUE, show_column_names = TRUE,  # Show row and column names
        rect_gp = gpar(col = "black"),  # Border color of cells
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", mat_df[i, j]), x = x, y = y, gp = gpar(fontsize = 10))
        }  # Display values within cells
)



pi1_anderson_leige = calculatePi1(table1 = all[["Georges_enterocyte_qval"]][,req_cols],
                                    table2 = all[["Anderson_enterocyte_qval"]][,req_cols],
                                    qvalue_thresh = 0.05,
                                    feature_id = "phenotype_id")