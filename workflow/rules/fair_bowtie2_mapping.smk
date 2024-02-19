module fair_bowtie2_mapping:
    snakefile:
        github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="3.1.0")
    config:
        config


use rule * from fair_bowtie2_mapping