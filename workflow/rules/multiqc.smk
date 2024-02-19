rule fair_gatk_mutect_somatic_multiqc_report:
    input:
        unpack(get_fair_gatk_mutect_somatic_multiqc_report_input),
    output:
         report(
            "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_Somatic.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "VariantCalling",
                "organism": "{species}.{build}.{release}.{datatype}",
            },
        ),
        "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_Mapping_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    params:
        extra=lookup(dpath="params/multiqc", within=config),
        use_input_files_only=True,
    log:
        "logs/fair_bowtie2_mapping/multiqc_report/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/multiqc_report/{species}.{build}.{release}.{datatype}.tsv"
    wrapper:
        "v3.3.6/bio/multiqc"