## Filter ABC predictions for E-G pairs overlapping CRISPR E-G pairs

# save.image("filt_abc.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# load ABC predictions
abc <- fread(snakemake@input$abc)

# load CRISPR data
crispr <- fread(snakemake@input$crispr)

# create GenomicRanges objects for ABC and CRISPR E-G pairs
crispr_gr <- with(crispr, GRanges(seqnames = paste0(chrom, ":", measuredGeneSymbol),
                                  ranges = IRanges(chromStart, chromEnd)))

abc_gr <- with(abc, GRanges(seqnames = paste0(chr, ":", TargetGene), ranges = IRanges(start, end)))

# set same seqlevels for both GRanges objects to avoid warnings
seqlevels_all_pairs <- as.character(unique(c(seqnames(crispr_gr), seqnames(abc_gr))))
seqlevels(crispr_gr) <- seqlevels_all_pairs
seqlevels(abc_gr) <- seqlevels_all_pairs

# extend all CRISPR enhancers by 1kb on each size to be sure to capture all CRISPR relevant
# enhancers
# start(crispr_gr) <- start(crispr_gr) - 1000
# end(crispr_gr) <- end(crispr_gr) + 1000

# merge now overlapping enhancers
# crispr_gr <- reduce(crispr_gr)

# find overlaps between ABC and CRISPR E-G pairs
ovl <- findOverlaps(abc_gr, crispr_gr)

# filter ABC predictions for E-G pairs overlapping CRISPR E-G pairs
abc_filt <- abc[unique(queryHits(ovl)), ]

# save to output file
fwrite(abc_filt, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
