## Convert element-level feature table to feature table using ABC E-G pair universe

# required packages
suppressPackageStartupMessages({
  library(data.table)
})

# load enhancer activity feature table
features <- fread(snakemake@input$features)

# load abc predictions
abc <- fread(snakemake@input$abc)

# extract ABC E-G pairs and merge with feature table based on enhancer coordinates
eg_pairs <- abc[, c("chr", "start", "end", "name", "TargetGene")]
features <- merge(eg_pairs, features, by = c("chr", "start", "end"), all.x = TRUE)

# write new feature table to output file
fwrite(features, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
