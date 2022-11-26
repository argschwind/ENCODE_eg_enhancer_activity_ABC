## Combine different enhancer activity measurements into one table

# save.image("combine.rda")
# stop()

# required packages
library(data.table)

# files containing assay quantifications from run.neighborhoods.py
infiles <- unlist(snakemake@input)
names(infiles) <- basename(dirname(infiles))

# load all quantification files
quants <- lapply(infiles, FUN = fread)

# only retain relevant columns and rename quantification column to correct assay
quants <- lapply(names(quants), FUN = function(assay){
  cols <- c("chr", "start", "end", "class", "name", "normalized_h3K27ac", "activity_base")
  output <- quants[[assay]][, ..cols]
  colnames(output)[6] <- paste0("normalized.", tolower(assay), "_quant")
  colnames(output)[7] <- paste0("activity_base.", tolower(assay), "_quant")
  return(output)
})

# merge all features into one table
quants <- Reduce(
  function(df1, df2) {
    merge(df1, df2, by = c("chr", "start", "end", "class", "name"), all = TRUE)
  }, quants)

# save to output file
fwrite(quants, file = snakemake@output[[1]], sep = "\t")
