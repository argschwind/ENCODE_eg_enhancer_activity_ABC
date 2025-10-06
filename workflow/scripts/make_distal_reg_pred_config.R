## Create pred_config file for CRISPR benchmarking

# save.image("distal_pred_config.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# load enhancer activity ABC models
pred_assayOnly <- fread(snakemake@input$pred_config_assayOnly)
pred_assayDHS <- fread(snakemake@input$pred_config_assayDHS)

# load distal regulation pred_config file containing config for other predictive models and retain
# columns relevant for CRISPR benchmarking
distal_reg_config <- fread(snakemake@input$pred_config)
distal_reg_config <- distal_reg_config[, 1:10]

# extract full ABC, 3D contact and ENCODE-E2Gextended from distal regulation pred_config
distal_reg_config <- distal_reg_config %>% 
  filter(pred_id %in% c("ABCdnase", "ABCfull", "ENCODE_rE2G", "ENCODE_rE2Gext", "HiC")) %>% 
  mutate(pred_id = replace(pred_id, pred_id == "HiC", "ABCfull"),
         pred_col = replace(pred_col, pred_col == "hic_contact", "hic_contact_pl_scaled_adj"),
         plot = TRUE)

# extract predictor score for DHS x H3K27ac (ABC) model from assayDHS pred_config table
dhs_pred_config <- filter(pred_assayDHS, grepl(pred_name_long, pattern = "H3K27ac"))

# combine with assayOnly models to create output
pred_config <- bind_rows(distal_reg_config, dhs_pred_config, pred_assayOnly)
  
# save pred_config to output file
write_tsv(pred_config, file = snakemake@output[[1]])
