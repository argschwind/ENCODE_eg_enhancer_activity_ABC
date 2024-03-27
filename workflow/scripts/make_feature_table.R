## Combine ENCODE enhancer activity assay quantifications to build feature table for all assays

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load ENCODE enhancer activity assay metadata
meta <- fread(snakemake@input$meta, na.strings = "")

# add column uniquely identifying the assay in each file
meta <- meta %>% 
  mutate(assay_uid = if_else(is.na(`Experiment target`), true = Assay,
                             false = `Experiment target`)) %>% 
  mutate(assay_uid = paste0(sub("-human", "", assay_uid), ".", `File accession`))

# files containing activity quantifications for each assay
assay_files <- unlist(snakemake@input[snakemake@input != snakemake@input$meta])
names(assay_files) <- basename(dirname(assay_files))

# load all assay files and retain relevant columns and combine into one table in long format
assays <- lapply(assay_files, FUN = fread, select = c("chr", "start", "end", "normalized_h3K27ac"))
assays <- rbindlist(assays, idcol = "File accession")

# add unique identifier for assays to table
assays <- assays %>% 
  left_join(select(meta, `File accession`, assay_uid), by = "File accession") %>% 
  select(-`File accession`)

# convert to wide format, one column per assay
assays <- dcast(assays, chr + start + end ~ assay_uid, value.var = "normalized_h3K27ac")

# write to output file
fwrite(assays, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
