module gatk_mutect2_calling:
    meta_wrapper:
        "v3.3.6/meta/bio/gatk_mutect2_calling"
    config:
        config



use rule picard_replace_read_groups from gatk_mutect2_calling as fair_gatk_mutect_somatic_picard_replace_read_groups with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        temp(
            "tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="--RGPL illumina --RGLB WES_Germline_{species}.{build}.{release}.{datatype} --RGPU {sample} --RGSM {sample}",


use rule sambamba_index_picard_bam from gatk_mutect2_calling as fair_gatk_mutect_somatic_sambamba_index_picard_bam with:
    input:
        "tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/sambamba_index_picard_bam/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/sambamba_index_picard_bam/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.tsv"


use rule mutect2_call from gatk_mutect2_calling as fair_gatk_mutect_somatic_gatk_mutect2 with:
    input:
        unpack(get_fair_gatk_mutect_somatic_gatk_mutect2_input),
    output:
        vcf=temp("tmp/fair_gatk_mutect_somatic/gatk_mutect2/{species}.{build}.{release}.{datatype}/{sample}.vcf"),
        f1r2=temp(
            "tmp/fair_gatk_mutect_somatic/gatk_mutect2/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/gatk_mutect2/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/gatk_mutect2/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        use_parallelgc=True,
        use_omp=True,
        extra=lookup(dpath="params/gatk/mutect2", within=config),


use rule gatk_get_pileup_summaries from gatk_mutect2_calling as fair_gatk_mutect_somatic_gatk_get_pileup_summaries with:
    input:
        unpack(get_fair_gatk_mutect_somatic_gatk_get_pileup_summaries_input),
    output:
        temp(
            "tmp/fair_gatk_mutect_somatic/gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("gatk", {}).get("getpileupsummaries", ""),


use rule gatk_calculate_contamination from gatk_mutect2_calling as fair_gatk_mutect_somatic_gatk_calculate_contamination with:
    input:
        "tmp/fair_gatk_mutect_somatic/gatk_get_pileup_summariess/{species}.{build}.{release}.{datatype}/{sample}.table",
    output:
        temp(
            "tmp/fair_gatk_mutect_somatic/gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/gatk/calculatecontamination", within=config),


use rule gatk_learn_read_orientation_model from gatk_mutect2_calling as fair_gatk_mutect_somatic_gatk_learn_read_orientation_model with:
    input:
        f1r2="tmp/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz",
    output:
        temp(
            "tmp/fair_gatk_mutect_somatic/gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.tar.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/gatk/learnreadorientationmodel", within=config),


use rule filter_mutect_calls from gatk_mutect2_calling as fair_gatk_mutect_somatic_gatk_filter_mutect_calls with:
    input:
        unpack(get_fair_gatk_mutect_somatic_gatk_filter_mutect_calls_input),
    output:
        vcf="results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz",
        vcf_tbi="results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_somatic/gatk_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_somatic/gatk_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/gatk/filtermutectcalls", within=config),