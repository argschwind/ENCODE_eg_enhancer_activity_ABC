## Rules to build ABC variations using different chromatin assays to estimate activity

# define python input functions
def get_default_dnase_bam(wildcards):
  dnase_id = config["default_dnase"][wildcards.cell_type]["bam"]
  dnase_file = config["scratch"] + "/enhancer_activity/bam/" + wildcards.cell_type + "/DNase-seq/" + dnase_id + "/" + dnase_id + ".sorted.bam"
  return(dnase_file)

def get_default_dnase_bam_bai(wildcards):
  dnase_id = config["default_dnase"][wildcards.cell_type]["bam"]
  dnase_file = config["scratch"] + "/enhancer_activity/bam/" + wildcards.cell_type + "/DNase-seq/" + dnase_id + "/" + dnase_id + ".sorted.bam.bai"
  return(dnase_file)

def get_default_dnase_bigwig(wildcards):
  dnase_id = config["default_dnase"][wildcards.cell_type]["bigWig"]
  dnase_file = config["scratch"] + "/enhancer_activity/bigWig/" + wildcards.cell_type + "/DNase-seq/" + dnase_id + "/" + dnase_id + ".bigWig"
  return(dnase_file)

# Prepare input files ------------------------------------------------------------------------------

# create chromosome sizes bed file
rule create_chrom_sizes_bed_file:
	input: "workflow/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv"
	output: "resources/GRCh38_EBV.no_alt.chrom.sizes.bed"
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input} > {output}
		"""

# download specified CRISPR data
rule download_crispr:
  output: "resources/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
  params:
    url = config["crispr_url"]
  conda: "../envs/enhancer_activity.yml"
  shell:
    "wget -O {output} {params.url}"

# filter ABC predictions for E-G pairs overlapping CRISPR E-G pairs
rule filter_abc_crispr:
  input:
    abc = lambda wildcards: config["abc_predictions"][wildcards.cell_type],
    crispr = config["crispr_data"]
  output: "resources/{cell_type}_ABC_predictions.CRISPR_only.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "8G",
  script:
    "../scripts/filter_abc_crispr.R"
  
    
# Download data from ENCODE ------------------------------------------------------------------------

# download chromatin assay bam of bigWig files
rule download_assay_data:
  output: temp(config["scratch"] + "/enhancer_activity/{ext}/{cell_type}/{assay}/{file}/{file}.{ext}")
  params:
    url = lambda wildcards: config["assays"][wildcards.ext][wildcards.cell_type][wildcards.assay][wildcards.file]
  wildcard_constraints:
    ext="bam|bigWig"
  conda: "../envs/enhancer_activity.yml"
  shell:
    "wget -O {output} {params.url}"
    
# sort and index bam files
rule sort_bam:
  input: config["scratch"] + "/enhancer_activity/bam/{cell_type}/{assay}/{file}/{file}.bam"
  output:
    bam = temp(config["scratch"] + "/enhancer_activity/bam/{cell_type}/{assay}/{file}/{file}.sorted.bam"),
    bai = temp(config["scratch"] + "/enhancer_activity/bam/{cell_type}/{assay}/{file}/{file}.sorted.bam.bai")
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  shell:
    "samtools sort -o {output.bam} {input}; samtools index {output.bam}"
    

# Compute ABC scores for different assays ----------------------------------------------------------

# quantify assays in candidate regions using bam files as input
rule run_neighborhoods_bam:
  input:
    candidate_regions = lambda wildcards: config["abc_candidate_regions"][wildcards.cell_type],
    access_bam = get_default_dnase_bam,
    access_bai = get_default_dnase_bam_bai,
    activity_bam = config["scratch"] + "/enhancer_activity/bam/{cell_type}/{assay}/{file}/{file}.sorted.bam",
    activity_bai = config["scratch"] + "/enhancer_activity/bam/{cell_type}/{assay}/{file}/{file}.sorted.bam.bai",
    genes = config["genes"],
    chrs = "workflow/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv",
    chrs_bed = "resources/GRCh38_EBV.no_alt.chrom.sizes.bed",
    ubiq_genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenes.txt"
  output: config["scratch"] + "/enhancer_activity/ABC/{cell_type}/bam/{assay}/{file}/EnhancerList.txt"
  params:
    outdir = config["scratch"] + "/enhancer_activity/ABC/{cell_type}/bam/{assay}/{file}"
  conda: "../ABC-Enhancer-Gene-Prediction/workflow/envs/abcenv.yml"
  resources:
    mem = "32G",
    runtime = "2h"
  shell:
    "python workflow/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py "
    "--candidate_enhancer_regions {input.candidate_regions} "
    "--DHS {input.access_bam} "
    "--default_accessibility_feature DHS "
    "--H3K27ac {input.activity_bam} "
    "--outdir {params.outdir} "
    "--chrom_sizes {input.chrs} "
    "--chrom_sizes_bed {input.chrs_bed} "
    "--genes {input.genes} "
    "--ubiquitously_expressed_genes {input.ubiq_genes} "
    
# quantify assays in candidate regions using bigWig files as input
rule run_neighborhoods_bigwig:
  input:
    candidate_regions = lambda wildcards: config["abc_candidate_regions"][wildcards.cell_type],
    access_bigwig = get_default_dnase_bigwig,
    activity_bigwig = config["scratch"] + "/enhancer_activity/bigWig/{cell_type}/{assay}/{file}/{file}.bigWig",
    genes = config["genes"],
    chrs = "workflow/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv",
    chrs_bed = "resources/GRCh38_EBV.no_alt.chrom.sizes.bed",
    ubiq_genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenes.txt"
  output: config["scratch"] + "/enhancer_activity/ABC/{cell_type}/bigWig/{assay}/{file}/EnhancerList.txt"
  params:
    outdir = config["scratch"] + "/enhancer_activity/ABC/{cell_type}/bigWig/{assay}/{file}"
  conda: "../ABC-Enhancer-Gene-Prediction/workflow/envs/abcenv.yml"
  resources:
    mem = "32G",
    runtime = "2h"
  shell:
    "python workflow/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py "
    "--candidate_enhancer_regions {input.candidate_regions} "
    "--DHS {input.access_bigwig} "
    "--default_accessibility_feature DHS "
    "--H3K27ac {input.activity_bigwig} "
    "--outdir {params.outdir} "
    "--chrom_sizes {input.chrs} "
    "--chrom_sizes_bed {input.chrs_bed} "
    "--genes {input.genes} "
    "--ubiquitously_expressed_genes {input.ubiq_genes} "

# compute ABC scores using different enhancer activity assays
rule compute_abc_scores:
  input:
    quants = config["scratch"] + "/enhancer_activity/ABC/{cell_type}/{file_type}/{assay}/{file}/EnhancerList.txt",
    abc = lambda wildcards: config["abc_predictions"][wildcards.cell_type],
    abc_crispr = "resources/{cell_type}_ABC_predictions.CRISPR_only.tsv.gz"
  output: 
    abc_scores_full = config["scratch"] + "/enhancer_activity/ABC/{cell_type}/{file_type}/{assay}/{file}/{type}_abc_scores_full.txt",
    abc_scores_crispr = config["scratch"] + "/enhancer_activity/ABC/{cell_type}/{file_type}/{assay}/{file}/{type}_abc_scores_crispr.txt"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/compute_abc.R"

# combine enhancer activity ABC scores into one table
rule combine_abc_scores:
  input:
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/K562/{{file_type}}/DNase-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type]["K562"]["DNase-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/K562/{{file_type}}/Histone_ChIP-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type]["K562"]["Histone_ChIP-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/K562/{{file_type}}/TF_ChIP-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type]["K562"]["TF_ChIP-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/K562/{{file_type}}/ATAC-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type]["K562"]["ATAC-seq"])
  output: temp("results/{file_type}/K562/{type}_abc_scores_{univ}.tsv.gz")
  conda: "../envs/enhancer_activity.yml"
  shell:
    "paste {input} | gzip > {output}"
    
# combine CRISPR ABC scores with E-G pairs to create predictions in ENCODE4 format
rule assemble_abc_predictions_crispr:
  input:
    abc = "resources/{cell_type}_ABC_predictions.CRISPR_only.tsv.gz",
    abc_scores = "results/{file_type}/{cell_type}/{type}_abc_scores_crispr.tsv.gz"
  output: "results/{file_type}/{cell_type}/{type}_abc_models_crispr.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/assemble_abc_predictions.R"

# combine genome-wide ABC scores with E-G pairs to create predictions in ENCODE4 format
rule assemble_abc_predictions_full:
  input:
    abc = lambda wildcards: config["abc_predictions"][wildcards.cell_type],
    abc_scores = "results/{file_type}/{cell_type}/{type}_abc_scores_full.tsv.gz"
  output: "results/{file_type}/{cell_type}/{type}_abc_models_full.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "48G"
  script:
    "../scripts/assemble_abc_predictions.R"    

# Create config files for CRISPR benchmarking ------------------------------------------------------
    
# create pred_config file for CRISPR benchmarking
rule make_pred_config:
  input:
    pred = "results/{file_type}/{cell_type}/{pred_type}_abc_models_crispr.tsv.gz",
    metadata = "resources/processed_encode_chromatin_metadata_{file_type}.tsv.gz",
  output: temp("results/{file_type}/{cell_type}/{pred_type}_pred_config.tsv")
  wildcard_constraints:
    pred_type = "assayOnly|assayDHS"
  params:
    default_dnase = lambda wildcards: config["default_dnase"][wildcards.cell_type][wildcards.file_type]
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/make_pred_config.R"
    
# create pred_config file for CRISPR benchmarking including other distal regulation predictors
rule make_distal_reg_pred_config:
  input:
    pred_config_assayOnly = "results/{file_type}/K562/assayOnly_pred_config.tsv",
    pred_config_assayDHS = "results/{file_type}/K562/assayDHS_pred_config.tsv",
    pred_config = config["distal_reg_pred_config"]
  output: "results/{file_type}/K562/EnhActABC_distal_reg_pred_config_{file_type}.tsv"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/make_distal_reg_pred_config.R"
