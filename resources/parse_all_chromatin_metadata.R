## Query the ENCODE portal for K562 chromatin data, download results as metadata table and filter
## to obtain one bam file per filtered experiment. Save processed metadata to table and create .yaml
## file containing download URLs for each file to load into snakemake config object.

# required packages
library(data.table)
library(tidyverse)
library(here)
library(yaml)

# Download ENCODE metadata -------------------------------------------------------------------------

# elements of metadata search URL
base = "https://www.encodeproject.org/metadata/?"
dataset_type = "type=Experiment&control_type%21=%2A&status=released&perturbed=false"
biosample = "biosample_ontology.term_name=K562"
assays = c("TF+ChIP-seq", "Histone+ChIP-seq", "ATAC-seq", "DNase-seq")
genome = "assembly=GRCh38"
file_type = "files.file_type=bam"
released = "files.analyses.status=released"

# assemble URL
assays <- paste0("assay_title=", paste(assays, collapse = "&assay_title="))
url <- paste0(base, paste(dataset_type, biosample, assays, genome, file_type, released, sep = "&"))

# query and download metadata as table
metadata <- fread(url)

# save complete metadata to file
metadata_file <- here("resources/all_k562_chromatin_metadata.tsv.gz")
fwrite(metadata, file = metadata_file, sep = "\t", quote = FALSE, na = "NA")

# General filters ----------------------------------------------------------------------------------

# only keep GRCh38 files
metadata <- filter(metadata, `File assembly` == "GRCh38")

# process 'DNase-seq' ------------------------------------------------------------------------------

dnase_seq <- metadata %>%
  filter(Assay == "DNase-seq") %>%
  filter(`Output type` == "unfiltered alignments" & `File type` == "bam") %>%
  filter(`Technical replicate(s)` == "2_1") %>%
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# if there are two experiments per target, pick the newer experiment or all if all have same date
dnase_seq <- dnase_seq %>%
  group_by(`Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

# process 'Histone ChIP-seq' -----------------------------------------------------------------------

# get all histone ChIP-seq files for technical replicate 2_1 for each experiment
histone_chipseq <- metadata %>%
  filter(Assay == "Histone ChIP-seq") %>%
  filter(`Output type` == "unfiltered alignments" & `File type` == "bam") %>%
  filter(`Technical replicate(s)` == "2_1") %>%
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# if there are two experiments per target, pick the newer experiment or all if all have same date
histone_chipseq <- histone_chipseq %>%
  group_by(`Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

# process 'TF ChIP-seq' ----------------------------------------------------------------------------

tf_chipseq <- metadata %>%
  filter(Assay == "TF ChIP-seq") %>%
  filter(`Output type` == "unfiltered alignments" & `File type` == "bam") %>%
  filter(`Technical replicate(s)` == "2_1") %>%
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# if there are two experiments per target, pick the newer experiment or all if all have same date
tf_chipseq <- tf_chipseq %>%
  group_by(`Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

# process 'ATAC-seq' -------------------------------------------------------------------------------

atac_seq <- metadata %>%
  filter(Assay == "ATAC-seq") %>%
  filter(`Output type` == "unfiltered alignments" & `File type` == "bam") %>%
  filter(`Technical replicate(s)` == "2_1") %>%
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# if there are two experiments per target, pick the newer experiment or all if all have same date
atac_seq <- atac_seq %>%
  group_by(`Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

# combine into one table ---------------------------------------------------------------------------

# combine metadata tables for individual assays
processed_metadata <- bind_rows(dnase_seq, histone_chipseq, tf_chipseq, atac_seq)

# save full processed metadata table to file
proc_metadata_file <- here("resources/processed_k562_chromatin_metadata.tsv.gz")
fwrite(processed_metadata, file = proc_metadata_file, sep = "\t", quote = FALSE, na = "NA")

# get file download URLs for each experiment in processed metadata
sample_urls <- processed_metadata %>%
  ungroup() %>%
  select(Assay, `File accession`, `File download URL`) %>% 
  mutate(Assay = gsub(" ", "_", Assay))  # replace spaces by underscores for simpler file paths

# convert to list of lists, containing download URLs for each file per assay
sample_urls <- sample_urls %>%
  as.data.table() %>%
  split(., by = "Assay", keep.by = FALSE) %>%
  lapply(., FUN = function(x) as.list(deframe(x)) )

# create yaml file with download urls per experiment to add to snakemake workflow
sample_urls <- list(assays = sample_urls)
write_yaml(sample_urls, file = here("config/encode4_k562_chrom_data.yml"))
