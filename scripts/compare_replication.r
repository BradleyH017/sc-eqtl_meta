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
library(circlize) 
library(grid)
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
    if("pval_perm" %in% colnames(all[[i]]) & !("p_beta" %in% colnames(all[[i]]))){
        colnames(all[[i]])[which(colnames(all[[i]]) == "pval_perm")] = "p_beta" # CHECK THIS (APPLIES TO FRANKE DATASET OF SIG QVALS)
    }
}

################## 1. Look at the overlap of genes tested ###################
get_sig_hits = function(x, thresh){
    return(unique(x[x$qvalue < thresh,]$phenotype_id))
}

# Setting no threshold here
list_data <- list(Georges = get_sig_hits(all[["Georges_enterocyte_qval"]], 1), 
                    Georges_tier2 = get_sig_hits(all[["Georges_enterocyte_tier_2_qval"]], 1), 
                    Anderson = get_sig_hits(all[["Anderson_enterocyte_qval"]], 1), Franke = unique(all$Franke_enterocyte_sumstat_qval$phenotype_id))

unique_ids <- unique(unlist(list_data))
df <- sapply(list_data, function(set) unique_ids %in% set)
df <- as.data.frame(df)
rownames(df) <- unique_ids
df_long <- df %>%
  mutate(ID = rownames(df)) %>%
  pivot_longer(-ID, names_to = "Set", values_to = "Present")

# Create the upset plot using ComplexUpset
up = upset(df, colnames(df), encode_sets=T)
ggsave(filename = paste0(pathOut, "/upset_qval_tested_genes.png"), plot = up, width = 10, height = 6)


################## 2. Look at the overlap of significant eGenes, irrespective of direction ###################
# repeat above but with a threshold
list_data <- list(Georges = get_sig_hits(all[["Georges_enterocyte_qval"]], 0.05), 
                    Georges_tier2 = get_sig_hits(all[["Georges_enterocyte_tier_2_qval"]], 0.05), 
                    Anderson = get_sig_hits(all[["Anderson_enterocyte_qval"]], 0.05), Franke = unique(all$Franke_enterocyte_qval$phenotype_id))

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

all_use <- keep(all, ~ contains_all_req_cols(.x, req_cols))


# generate a pairwise matrix for replication in one dataset vs another
plot_names = c("Georges", "Georges_tier2", "Anderson", "Franke")
pigrid = matrix(nrow=length(all_use), ncol = length(all_use))
rownames(pigrid) = plot_names
colnames(pigrid) = plot_names

for ( r in 1:length(all_use) ){
    for ( c in 1:length(all_use) ){
        if(c != r){
            try({
                # Ref dataset is the row
                print(paste0("Ref ", names(all_use)[r], " vs query ", names(all_use)[c]))
                pi1 = calculatePi1(table1 = all_use[[r]][,req_cols],
                                        table2 = all_use[[c]][,req_cols],
                                        qvalue_thresh = 0.1,
                                        feature_id = "phenotype_id")
                # Add
                pigrid[r,c] = pi1
            })
        }
    }
}

# Plot this:
hm <- Heatmap(pigrid,
        name = "pi1 (1-pi0)",  # Name of the matrix
        cluster_rows=F,
        cluster_columns=F,
        col = colorRamp2(c(0, 1), c("white", "red")),  # Color scheme (adjust as needed)
        show_row_names = TRUE, show_column_names = TRUE,  # Show row and column names
        rect_gp = gpar(col = "black"),  # Border color of cells
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", pigrid[i, j]), x = x, y = y, gp = gpar(fontsize = 10))
        },  # Display values within cells
        column_names_rot = 45,
        row_title = "Reference",
        column_title = "Query",  # X-axis label
        heatmap_legend_param = list(title = "pi1 (1-pi0)", 
                                    title_position = "topcenter", 
                                    legend_direction = "horizontal")
)

ppi=300
png(paste0(pathOut, "/heatmap_pi1.png"), width = 9*ppi, height = 8*ppi, res=ppi)
draw(hm, padding = unit(c(3, 3, 2, 2), "cm"), 
        heatmap_legend_side = "right",              # Position legend on the right
        annotation_legend_side = "right" )
dev.off()

################## 4. Looking at the correlation of the direction of affects ###################
# For each dataset, select the top hit per phenotype (i.e lowest qval)
# This is already done in the anderson and Georges datasets
process_dataframe <- function(df) {
    df$qvalue = as.numeric(df$qvalue)
    df %>%
        group_by(phenotype_id) %>%
        arrange(qvalue) %>%
        slice_head(n = 1) %>% 
        ungroup() 
}

# Apply the function to each dataframe in the list
all_use_dir <- map(all_use, process_dataframe)

# Look at the pairwise direction of affects across datasets
beta_plots = list()
for ( r in 1:length(all_use_dir) ){
    for ( c in 1:length(all_use_dir) ){
        if(c > r){
            # Get hits
            print(paste0(plot_names[r], " vs ", plot_names[c]))
            common = intersect(all_use_dir[[r]]$phenotype_id, all_use_dir[[c]]$phenotype_id)
            tempr = all_use_dir[[r]][all_use_dir[[r]]$phenotype_id %in% common,]
            tempr = tempr[order(tempr$phenotype_id),]
            tempc = all_use_dir[[c]][all_use_dir[[c]]$phenotype_id %in% common,]
            tempc = tempc[order(tempc$phenotype_id),]
            sig1 = tempr[tempr$qvalue < 0.05,]$phenotype_id
            sig2 = tempc[tempc$qvalue < 0.05,]$phenotype_id
            sigboth = intersect(sig1, sig2)
            sig1 = setdiff(sig1, sigboth)
            sig2 = setdiff(sig2, sigboth)
            merged_data <- inner_join(tempr, tempc, by = "phenotype_id", suffix = c("_r", "_c"))
            merged_data$sig_group <- ifelse(merged_data$phenotype_id %in% sigboth, "sigboth",
                                ifelse(merged_data$phenotype_id %in% sig1, "sig1",
                                        ifelse(merged_data$phenotype_id %in% sig2, "sig2", "none")))

            # Get prportion per quadrant
            quad1_prop <- sum(merged_data$slope_r > 0 & merged_data$slope_c > 0) / length(common)
            quad2_prop <- sum(merged_data$slope_r > 0 & merged_data$slope_c < 0) / length(common)
            quad3_prop <- sum(merged_data$slope_r < 0 & merged_data$slope_c > 0) / length(common)
            quad4_prop <- sum(merged_data$slope_r < 0 & merged_data$slope_c < 0) / length(common)


            # plot
            p = ggplot(merged_data, aes(x = slope_r, y = slope_c, color = sig_group)) +
                        geom_point(size = 2.5, alpha = 0.7) +
                        scale_color_manual(values = c("sig1" = "orange", "sig2" = "navy", "sigboth" = "darkgreen", "none" = "gray"),
                                            labels = c("sig1" = "x only", "sig2" = "y only", "sigboth" = "both", "none" = "Neither")) +
                        labs(x = "tempr slope", y = "tempc slope", color = "Signature Group") +
                        theme_minimal() +
                        theme(panel.background = element_rect(fill = "white", color = "black"),
                        plot.background = element_rect(fill = "white", color = "black"),
                        plot.margin = margin(20, 20, 20, 20),
                        panel.border = element_rect(color = "black", fill = NA, size = 1)) +  
                        xlab(paste0("Effect size:", plot_names[r])) + 
                        ylab(paste0("Effect size:", plot_names[c])) + 
                        xlim(c(-1.75,1.75)) + 
                        ylim(c(-1.75,1.75)) +  
                        ggtitle(paste0(plot_names[r], " vs ", plot_names[c], ": qvalue < 0.05")) + 
                        geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
                        geom_vline(xintercept = 0, color = "black", linetype = "solid") +
                        theme(panel.border = element_rect(color = "black", fill = NA, size = 2)) + 
                        geom_text(aes(x = Inf, y = Inf, label = scales::percent(quad1_prop), color = "black"), hjust = 2.0, vjust = 2.0) +
                        geom_text(aes(x = Inf, y = -Inf, label = scales::percent(quad2_prop), color = "black"), hjust = 2.0, vjust = -1) +
                        geom_text(aes(x = -Inf, y = Inf, label = scales::percent(quad3_prop), color = "black"), hjust = -1, vjust = 2.0) +
                        geom_text(aes(x = -Inf, y = -Inf, label = scales::percent(quad4_prop), color = "black"), hjust = -1, vjust = -1)


            beta_plots[[length(beta_plots) + 1]] = p
            ggsave(filename = paste0(pathOut, "/", plot_names[r], "_vs_", plot_names[c], "_effect_sizes_qval.png"), plot = p, width = 8, height = 6, dpi=300)
        }
    }
}

