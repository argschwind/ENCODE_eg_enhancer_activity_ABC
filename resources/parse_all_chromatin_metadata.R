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
released = "files.analyses.status=released"

# assemble URLs for bam, fastq and bigWig files
assays <- paste0("assay_title=", paste(assays, collapse = "&assay_title="))
bam_url <- paste0(base, paste(dataset_type, biosample, assays, genome, "files.file_type=bam", released, sep = "&"))
fastq_url <- paste0(base, paste(dataset_type, biosample, assays, genome, "files.file_type=fastq", sep = "&"))
bigwig_url <- paste0(base, paste(dataset_type, biosample, assays, genome, "files.file_type=bigWig", released, sep = "&"))

# query and download metadata as table
bam_metadata <- fread(bam_url)
fastq_metadata <- fread(fastq_url)
bigwig_metadata <- fread(bigwig_url)

# Process metadata  --------------------------------------------------------------------------------

# only keep ENCODE4 GRCh38 bam files for one technical replicate
bam_metadata <- bam_metadata %>% 
  filter(`File assembly` == "GRCh38") %>% 
  filter(`Technical replicate(s)` == "2_1") %>% 
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# only keep ENCODE4 GRCh38 coverage bigWig files for one technical replicate
bigwig_metadata <- bigwig_metadata %>% 
  filter(`File assembly` == "GRCh38") %>% 
  filter(`Output type` %in% c("fold change over control", "read-depth normalized signal")) %>% 
  filter(`Technical replicate(s)` == "2_1") %>% 
  filter(grepl(`File analysis title`, pattern = "^ENCODE4"))

# get Run type information (single-end vs paired-end reads) for each experiment from fastq metadata
run_type <- fastq_metadata %>% 
  filter(`Technical replicate(s)` == "2_1") %>% 
  select(`Experiment accession`, `Run type`) %>% 
  distinct()

# for one specific experiment, remove wrong type TODO: FIX THIS USING 'DERIVED FROM'
run_type <- run_type %>% 
  filter(!(`Experiment accession` == "ENCSR668LDD" & `Run type` == "paired-ended"))

# add Run type to bam metadata files
bam_metadata <- bam_metadata %>% 
  select(-`Run type`) %>% 
  left_join(run_type, by = "Experiment accession")

# add Run type to bigWig metadata files
bigwig_metadata <- bigwig_metadata %>% 
  select(-`Run type`) %>% 
  left_join(run_type, by = "Experiment accession")

# select filtered or unfiltered bam files based on whether data has single-end or paired-end reads
bam_metadata <- bam_metadata %>% 
  filter( (`Output type` == "alignments" & `Run type` == "paired-ended") | 
            (`Output type` == "unfiltered alignments" & `Run type` == "single-ended") )

# if there are two experiments per target, pick the newer experiment or all if all have same date
bam_metadata <- bam_metadata %>%
  group_by(Assay, `Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

bigwig_metadata <- bigwig_metadata %>%
  group_by(Assay, `Experiment target`) %>%
  slice_max(`Experiment date released`, n = 1)

# Save processed metadata to output files ----------------------------------------------------------

# save full processed metadata tables to file
bam_metadata_file <- here("resources/processed_k562_chromatin_metadata_bam.tsv.gz")
fwrite(bam_metadata, file = bam_metadata_file, sep = "\t", quote = FALSE, na = "NA")

bigwig_metadata_file <- here("resources/processed_k562_chromatin_metadata_bigWig.tsv.gz")
fwrite(bigwig_metadata, file = bigwig_metadata_file, sep = "\t", quote = FALSE, na = "NA")

# Create yaml format files containing download URLs for each file ----------------------------------

# get file download URLs for each experiment in processed bam metadata
bam_sample_urls <- bam_metadata %>%
  ungroup() %>%
  select(Assay, `File accession`, `File download URL`) %>% 
  mutate(Assay = gsub(" ", "_", Assay))  # replace spaces by underscores for simpler file paths

# get file download URLs for each experiment in processed bigWig metadata
bigwig_sample_urls <- bigwig_metadata %>%
  ungroup() %>%
  select(Assay, `File accession`, `File download URL`) %>% 
  mutate(Assay = gsub(" ", "_", Assay))  # replace spaces by underscores for simpler file paths

# convert to list of lists, containing download URLs for each file per assay
bam_sample_urls <- bam_sample_urls %>%
  as.data.table() %>%
  split(., by = "Assay", keep.by = FALSE) %>%
  lapply(., FUN = function(x) as.list(deframe(x)) )

# convert to list of lists, containing download URLs for each file per assay
bigwig_sample_urls <- bigwig_sample_urls %>%
  as.data.table() %>%
  split(., by = "Assay", keep.by = FALSE) %>%
  lapply(., FUN = function(x) as.list(deframe(x)) )

# combine into one list of lists
sample_urls <- list(assays = list(bam = bam_sample_urls, bigWig = bigwig_sample_urls))

# create yaml files with download urls per experiment to add to snakemake workflow
write_yaml(sample_urls, file = here("config/encode4_k562_chrom_data.yml"))
