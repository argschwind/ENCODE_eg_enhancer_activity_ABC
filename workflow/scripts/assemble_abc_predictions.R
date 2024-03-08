## Combine ABC scores with E-G pairs to create predictions in ENCODE4 format

# save.image("combine.rda")
# stop()

# required packages
suppressMessages({
  library(data.table)
  library(dplyr)
})

# load ABC file and retain only relevant columns
abc_cols <- c("chr", "start", "end", "name", "class", "TargetGene", "distance", "CellType")
abc <- fread(snakemake@input$abc, select = abc_cols)

# load enhancer activity ABC scores
abc_scores <- fread(snakemake@input$abc_scores)

# make sure that ABC table is ordered by name and TargetGene
abc <- arrange(abc, name, TargetGene)

# add enhancer activity ABC scores to ABC E-G pairs to create output
output <- cbind(abc, abc_scores)

# reformat to ENCODE4 predictions formant for output
output <- output %>% 
  mutate(TargetGeneEnsemblID = NA_character_, TargetGeneTSS = NA_integer_) %>% 
  select(chr, start, end, name, class, TargetGene, TargetGeneEnsemblID, TargetGeneTSS, CellType,
         ends_with(".Score"), DistanceToTSS = distance)

# save to output file
fwrite(output, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
