## Create CRISPR benchamrking pred_config file for one model and file type combination

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# load enhancer activity ABC models
pred <- fread(snakemake@input$pred)
  
# load metadata for ENCODE chromatin data
metadata <- fread(snakemake@input$metadata)

# filter metadata for input cell type
metadata <- filter(metadata, `Biosample term name` == snakemake@wildcards$cell_type)

# extract all predictor score columns and create table
pred_config <- grep(colnames(pred), pattern = ".Score", value = TRUE) %>%
  tibble(pred_col = .)

# create table containing meaningful names for each pred_col
pred_names <- metadata %>%
  select(`File accession`, `Assay`, `Experiment target`) %>%
  mutate(`Experiment target` = sub("-human", "", `Experiment target`)) %>% 
  mutate(pred_name_long = if_else(`Experiment target` == "",
                                  true = paste0(Assay, " (", `File accession`, ")"),
                                  false = paste0(`Experiment target`, " (", `File accession`, ")"))) %>%
  mutate(pred_col = paste0(`File accession`, ".Score")) %>%
  select(pred_col, pred_name_long)

# indicate in name if activity was computed by taking the geometric mean with DNase-seq
if (snakemake@wildcards$pred_type == "assayDHS") {
  dhs <- snakemake@params$default_dnase
  pred_names$pred_name_long <- paste0("DHS (", dhs, ") x ", pred_names$pred_name_long)
}
  
# add pretty names to pred config
pred_config <- left_join(pred_config, pred_names, by = "pred_col")

# add other required columns to create full pred_config file
pred_config <- pred_config %>% 
  mutate(pred_id = if_else(snakemake@wildcards$pred_type == "assayOnly",
                           true = "EnhActAssayOnly",
                           false = "EnhActAssayDHS")) %>% 
  mutate(boolean = FALSE, alpha = NA, aggregate_function = "sum", fill_value = 0,
         inverse_predictor = FALSE, color = "black", plot = FALSE) %>% 
  select(pred_id, pred_col, boolean, alpha, aggregate_function, fill_value, inverse_predictor, 
         pred_name_long, color, plot)

# save pred_config to output file
write_tsv(pred_config, file = snakemake@output[[1]])
