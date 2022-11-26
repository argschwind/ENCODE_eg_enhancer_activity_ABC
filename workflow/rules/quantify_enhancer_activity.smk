## Rules to build ABC variations using different chromatin assays to estimate activity

# download chromatin data in bam files
rule download_bam:
  output: temp("resources/assays/{sample}/{assay}_reads.bam")
  params:
    url = lambda wildcards: config["chrom_assays"][wildcards.sample][wildcards.assay]
  conda: "../envs/enhancer_activity.yml"
  shell:
    "wget -O {output} {params.url}"
    
# sorte and index bam files
rule sort_bam:
  input: "resources/assays/{sample}/{assay}_reads.bam"
  output:
    bam = temp("resources/assays/{sample}/{assay}_reads.sorted.bam"),
    bai = temp("resources/assays/{sample}/{assay}_reads.sorted.bam.bai")
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  shell:
    "samtools sort -o {output.bam} {input}; samtools index {output.bam}"
    
# make sure that gene annotations are sorted and unique
rule unique_gene_annot:
  input: 
    genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.bed",
    chrs = "resources/hg38.chrom.sizes"
  output: "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.sorted.unique.bed"
  conda: "../ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/abcenv.yml"
  shell:
    "bedtools sort -faidx {input.chrs} -i {input.genes} | uniq > {output}"

# quantify assays in candidate regions
rule run_neighborhoods:
  input:
    candidate_regions = lambda wildcards: config["abc_candidate_regions"][wildcards.sample],
    access_bam = "resources/assays/{sample}/DNase_reads.sorted.bam",
    access_bai = "resources/assays/{sample}/DNase_reads.sorted.bam.bai",
    activity_bam = "resources/assays/{sample}/{assay}_reads.sorted.bam",
    activity_bai = "resources/assays/{sample}/{assay}_reads.sorted.bam.bai",
    chrs = "resources/hg38.chrom.sizes",
    genes = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.sorted.unique.bed",
    ubiq_genes = "workflow/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt",
 #   qnorm = "workflow/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/EnhancersQNormRef.K562.txt"
  output: 
    enhancers = "results/{sample}/{assay}/EnhancerList.txt",
    genelist = "results/{sample}/{assay}/GeneList.txt"
  params:
    outdir = "results/{sample}/{assay}"
  conda: "../ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/abcenv.yml"
  resources:
    mem = "32G"
  shell:
    "python workflow/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/run.neighborhoods.py "
    "--candidate_enhancer_regions {input.candidate_regions} "
    "--DHS {input.access_bam} "
    "--H3K27ac {input.activity_bam} "
    "--outdir {params.outdir} "
    "--chrom_sizes {input.chrs} "
    "--genes {input.genes} "
    "--ubiquitously_expressed_genes {input.ubiq_genes} "
#    "--qnorm {input.qnorm}"

# combine quantifications for all assays per sample
rule combine_quantifications:
  input:
    lambda wildcards: expand("results/{{sample}}/{assay}/EnhancerList.txt",
      assay = config["chrom_assays"][wildcards.sample])
  output: "results/{sample}/activity_assay_quantifications.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/combine_quantifications.R"

# create feature table by merging enhancer quantifications with ABC table
rule create_feature_table:
  input: 
    quants = "results/{sample}/activity_assay_quantifications.tsv.gz",
    abc = lambda wildcards: config["abc_predictions"][wildcards.sample]
  output: "results/{sample}/enhancer_activity_feature_table.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  script:
    "../scripts/create_feature_table.R"
    
# compute ABC models usign different chromatin assays
rule chrom_abc_models:
  input: "results/{sample}/enhancer_activity_feature_table.tsv.gz"
  output: "results/{sample}/chrom_abc_models.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "48G"
  script:
    "../scripts/chrom_abc_models.R"
