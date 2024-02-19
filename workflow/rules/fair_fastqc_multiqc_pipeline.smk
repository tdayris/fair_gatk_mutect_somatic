module fair_fastqc_multiqc:
    snakefile:
        github("tdayris/fair_fastqc_multiqc", path="workflow/Snakefile", tag="2.0.2")
    config:
        config


use rule * from fair_fastqc_multiqc