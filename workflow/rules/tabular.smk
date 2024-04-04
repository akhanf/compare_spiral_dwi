rule import_hcp_demog:
    input:
        ndar_txt=lambda wildcards: config["in_demog_tabular"][wildcards.dataset],
        subjects_tsv="resources/dataset-{dataset}_subjects.tsv",
    output:
        tsv="resources/dataset-{dataset,HCP}_demographics.tsv",
    script:
        "../scripts/import_hcp_demog.py"
