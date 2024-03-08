## Compute ABC scores for a given enhancer activity assay. Computes ABC score using the assay alone
## to estimate enhancer activity and by using geometric mean together with DNase-seq

# save.image("compute_abc.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load regular ABC predictions and retain required columns
abc_cols <- c("name", "TargetGene", "hic_contact_pl_scaled_adj")
abc <- fread(snakemake@input$abc, select = abc_cols)

# load E-G pairs of ABC predictions overlapping CRISPR data
abc_crispr <- fread(snakemake@input$abc_crispr, select = c("name", "TargetGene"))

# columns to load from quantification files based on model 'type' wildcard
type <- snakemake@wildcards$type
if (type == "assayOnly") {
  quant_cols <- c("name", "normalized_h3K27ac")
} else if (type == "assayDHS") {
  quant_cols <- c("name", "activity_base")
} else {
  stop("Invalid 'type' wildcard. Must be either 'assayOnly' or 'assayDHS'", call. = FALSE)
}

# load relevant columns from quantification file
quant <- fread(snakemake@input$quants, select = quant_cols)
colnames(quant)[[2]] <- "activity"

# add quantifications to abc E-G pairs
abc <- merge(abc, quant, by = "name")

# compute abc numerators
abc$abc_numerator <- abc$activity * abc$hic_contact_pl_scaled_adj

# compute abc denominator TODO: replace by data.tables aggregate
abc <- abc %>% 
  group_by(TargetGene) %>% 
  mutate(abc_denominator = sum(abc_numerator, na.rm = TRUE))

# compute ABC score
abc$abc_score <- abc$abc_numerator / abc$abc_denominator

# make sure that ABC table is ordered by name and TargetGene remove grouping
abc <- abc %>% 
  ungroup() %>% 
  arrange(name, TargetGene)

# set new name column for ABC score which includes the chromatin file id
score_col <- sym(paste0(snakemake@wildcards$file, ".Score"))
abc <- rename(abc, !!score_col := abc_score)

# extract ABC scores for E-G pairs overlapping CRISPR E-G pairs
abc_crispr <- abc_crispr %>% 
  left_join(abc, by = c("name", "TargetGene")) %>% 
  arrange(name, TargetGene)

# save ABC scores for all E-G pairs to output file
abc %>% 
  select(all_of(score_col)) %>% 
  fwrite(file = snakemake@output$abc_scores_full, sep = ",", quote = FALSE, na = "NA")

# save ABC scores for CRISPR E-G pairs to output file
abc_crispr %>% 
  select(all_of(score_col)) %>% 
  fwrite(file = snakemake@output$abc_scores_crispr, sep = ",", quote = FALSE, na = "NA")
