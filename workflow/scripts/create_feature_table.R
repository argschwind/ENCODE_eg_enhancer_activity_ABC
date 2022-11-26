## Combine enhancer activity quantifications with abc predictions

# save.image("feature_table.rda")
# stop()

# required packages
library(data.table)

# load quantifications
quants <- fread(snakemake@input$quants)

# load abc predictions
abc <- fread(snakemake@input$abc)

# only retain relevant columns from quants
quants <- quants[, -c("chr", "start", "end", "class")]

# add quantifications to abc table
output <- merge(x = abc, y = quants, by = "name", all.x = TRUE)

# write output to file
fwrite(output, file = snakemake@output[[1]], sep = "\t")
