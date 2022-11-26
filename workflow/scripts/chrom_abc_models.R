## Compute ABC models using different chromatin assays to estimate enhancer activity

# save.image("abc.rda")
# stop()

# required packages
library(data.table)
library(tidyverse)

# load feature table
features <- fread(snakemake@input[[1]])

# get all columns containing chromatin features
chrom_feature_cols <- grep(colnames(features), pattern = "^.+_quant$", value = TRUE)

# convert to long format and reformat assay names
abc <- features %>% 
 select(name, TargetGene, normalized_dhs, powerlaw_contact, hic_contact_pl_scaled_adj,
        all_of(chrom_feature_cols)) %>%
  pivot_longer(cols = all_of(chrom_feature_cols), names_to = c("activity_type", "assay"),
               names_sep = "\\.", values_to = "activity") %>% 
  mutate(assay = sub("_quant", "", assay))

# compute abc numerators
abc <- abc %>% 
  mutate(abc_numerator = activity * hic_contact_pl_scaled_adj)

# compute abc denominator
abc <- abc %>% 
  group_by(TargetGene, activity_type, assay) %>% 
  mutate(abc_denominator = sum(abc_numerator))

# compute ABC score
abc <- abc %>% 
  mutate(abc_score = abc_numerator / abc_denominator)

# create score name based on activity type
abc <- abc %>% 
  mutate(
    score_name = if_else(activity_type == "activity_base",
                         true = paste0(assay, "_dhs.Score"),
                         false = paste0(assay, ".Score") )
  )
         
# reformat and convert to wide format
abc <- abc %>% 
  ungroup() %>% 
  select(name, TargetGene, score_name, abc_score) %>% 
  pivot_wider(names_from = score_name, values_from = abc_score)

# compute ABC scores using only one assay + DNase-seq ----------------------------------------------

# merge with other ABC columns from input to create output
output <- features %>% 
  select(-all_of(chrom_feature_cols)) %>% 
  left_join(abc, by = c("name", "TargetGene"))

# save to file
fwrite(output, file = snakemake@output[[1]], sep = "\t")
