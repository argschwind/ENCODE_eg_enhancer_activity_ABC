
# combine assay quantifications to create feature table for all ENCODE enhancer assays
rule make_feature_table:
  input:
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/{{cell_type}}/{{file_type}}/DNase-seq/{file}/EnhancerList.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type][wildcards.cell_type]["DNase-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/{{cell_type}}/{{file_type}}/Histone_ChIP-seq/{file}/EnhancerList.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type][wildcards.cell_type]["Histone_ChIP-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/{{cell_type}}/{{file_type}}/TF_ChIP-seq/{file}/EnhancerList.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type][wildcards.cell_type]["TF_ChIP-seq"]),
    lambda wildcards: expand("{scratch}/enhancer_activity/ABC/{{cell_type}}/{{file_type}}/ATAC-seq/{file}/EnhancerList.txt",
           scratch = config["scratch"], file = config["assays"][wildcards.file_type][wildcards.cell_type]["ATAC-seq"]),
    meta = "resources/processed_encode_chromatin_metadata_{file_type}.tsv.gz"
  output: temp("results/{file_type}/{cell_type}/enhancer_activity_features.tsv.gz")
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/make_feature_table.R"

# create feature table in ABC E-G universe
rule make_abc_feature_table:
  input: 
    features = "results/{file_type}/{cell_type}/enhancer_activity_features.tsv.gz",
    abc = lambda wildcards: config["abc_predictions"][wildcards.cell_type]
  output: "results/{file_type}/{cell_type}/EnhAct_features_ABC_egPairs.tsv.gz"
  conda: "../envs/enhancer_activity.yml"
  resources:
    mem = "48G"
  script:
    "../scripts/make_abc_feature_table.R"
