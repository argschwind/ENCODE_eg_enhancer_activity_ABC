## Rules to build ABC variations using different chromatin assays to estimate activity

# ADD: download chromosome sizes and extract regular chromosomes

# filter gene annotations for regular chromosomes only
rule filter_gene_annot:
  input:
    genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.bed",
    chrs = "resources/GRCh38_EBV.chrom.sizes.filt.bed"
  output: "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.regular_chroms.bed"
  conda: "../ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/abcenv.yml"
  shell:
    "bedtools intersect -a {input.genes} -b {input.chrs} -wa | "
    "bedtools sort -faidx {input.chrs} -i stdin > {output}"

# extract unique ABC candidate elements
rule abc_candidate_elements:
  input: 
    enh = lambda wildcards: config["abc_candidate_regions"][wildcards.celltype]
  output: "results/{celltype}/abc_candidate_elements.bed"
  conda: "../ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/abcenv.yml"
  shell:
    """awk 'BEGIN {{OFS = "\\t"}} NR>1 {{print $1, $2, $3}}' {input.enh} | uniq > {output}"""

# filter ABC predictions for E-G pairs overlapping CRISPR E-G pairs
rule filter_abc_crispr:
  input:
    abc = config["abc_predictions"]["K562"],
    crispr = config["crispr_data"]
  output: "resources/K562_ABC_predictions.CRISPR_only.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/filter_abc_crispr.R"

# download chromatin assay bigWig files
rule download_bigwig:
  output: temp(config["scratch"] + "/enhancer_activity/{assay}/{file}/{file}.bigWig")
  params:
    url = lambda wildcards: config["assays"][wildcards.assay][wildcards.file]
  conda: "../envs/enhancer_activity.yml"
  shell:
    "wget -O {output} {params.url}"
    
# # filter bigWig file for regular chromosomes and ensure that it's correctly sorted
# rule filter_bigwig:
#   input:
#     bigwig = config["scratch"] + "/enhancer_activity/{assay}/{file}/{file}.bigWig",
#     chrs = "resources/GRCh38_EBV.chrom.sizes.filt"
#   output: temp(config["scratch"] + "/enhancer_activity/{assay}/{file}/{file}.filt.bigWig")
#   conda: "../envs/enhancer_activity.yml"
#   script:
#     "../scripts/filter_bigwig.R"

# quantify assays in candidate regions
rule run_neighborhoods_bigwig:
  input:
    candidate_regions = "results/K562/abc_candidate_elements.bed",
    access_bigwig = config["scratch"] + "/enhancer_activity/DNase-seq/" + config["default_dnase"] + "/" + config["default_dnase"] + ".bigWig",
    activity_bigwig = config["scratch"] + "/enhancer_activity/{assay}/{file}/{file}.bigWig",
    chrs = "resources/GRCh38_EBV.chrom.sizes",
    chrs_bed = "resources/GRCh38_EBV.chrom.sizes.bed",
    genes = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.regular_chroms.bed",
    ubiq_genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt",
  output: config["scratch"] + "/enhancer_activity/{assay}/{file}/EnhancerList.txt"
  params:
    outdir = config["scratch"] + "/enhancer_activity/{assay}/{file}"
  conda: "../ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/abcenv.yml"
  resources:
    mem = "32G"
  shell:
    "python workflow/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/run.neighborhoods.py "
    "--candidate_enhancer_regions {input.candidate_regions} "
    "--DHS {input.access_bigwig} "
    "--H3K27ac {input.activity_bigwig} "
    "--outdir {params.outdir} "
    "--chrom_sizes {input.chrs} "
    "--genes {input.genes} "
    "--ubiquitously_expressed_genes {input.ubiq_genes} "

# compute ABC scores using different enhancer activity assays
rule compute_abc_scores:
  input:
    quants = config["scratch"] + "/enhancer_activity/{assay}/{file}/EnhancerList.txt",
    abc = config["abc_predictions"]["K562"],
    abc_crispr = "resources/K562_ABC_predictions.CRISPR_only.tsv.gz"
  output: 
    abc_scores_full = config["scratch"] + "/enhancer_activity/{assay}/{file}/{type}_abc_scores_full.txt",
    abc_scores_crispr = config["scratch"] + "/enhancer_activity/{assay}/{file}/{type}_abc_scores_crispr.txt"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/compute_abc.R"

# combine enhancer activity ABC scores into one table
rule combine_abc_scores:
  input:
    expand("{scratch}/enhancer_activity/DNase-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"]["DNase-seq"]),
    expand("{scratch}/enhancer_activity/Histone_ChIP-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"]["Histone_ChIP-seq"]),
    expand("{scratch}/enhancer_activity/TF_ChIP-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"]["TF_ChIP-seq"]),
    expand("{scratch}/enhancer_activity/ATAC-seq/{file}/{{type}}_abc_scores_{{univ}}.txt",
           scratch = config["scratch"], file = config["assays"]["ATAC-seq"])
  output: temp("results/K562/{type}_abc_scores_{univ}.tsv.gz")
  conda: "../envs/enhancer_activity.yml"
  shell:
    "paste {input} | gzip > {output}"
    
# combine ABC scores with E-G pairs to create predictions in ENCODE4 format
rule assemble_abc_predictions:
  input:
    abc = "resources/K562_ABC_predictions.CRISPR_only.tsv.gz",
    abc_scores = "results/K562/{type}_abc_scores_crispr.tsv.gz"
  output: "results/K562/{type}_abc_models_crispr.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/assemble_abc_predictions.R"
    
# create pred_config file for CRISPR benchmarkisng
rule make_pred_config:
  input:
    pred = "results/K562/{type}_abc_models_crispr.tsv.gz",
    metadata = "resources/processed_k562_chromatin_metadata_bigWig.tsv.gz"
  output: "results/K562/{type}_pred_config.tsv"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/make_pred_config.R"
    
# create pred_config file for CRISPR benchmarking including other distal regulation predictors
rule make_distal_reg_pred_config:
  input:
    pred_config_assayOnly = "results/K562/assayOnly_pred_config.tsv",
    pred_config_assayDHS = "results/K562/assayDHS_pred_config.tsv",
    pred_config = "/oak/stanford/groups/engreitz/Projects/Benchmarking/Predictors/benchmarking_pred_config.tsv"
  output: "results/K562/EnhActABC_distal_reg_pred_config.tsv"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/make_distal_reg_pred_config.R"
