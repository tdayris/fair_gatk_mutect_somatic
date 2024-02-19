import csv
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from typing import Any, NamedTuple

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.readline())
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = lookup(dpath="genomes", within=config)
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.readline())
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "transcripts"]
stream_list: list[str] = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


def get_fair_gatk_mutect_somatic_gatk_mutect2_input(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
    samples: pandas.DataFrame = samples
) -> dict[str, str]:
    """
    Return best input files list, according to
    Mutect2-snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)
    samples   (pandas.DataFrame)      : Describe samples

    Return (dict[str, str]):
    Dictionary of input files
    """
    mutect2_call_input: dict[str, str | List[str]] = {
        "fasta": None,
        "fasta_fai": None,
        "fasta_dict": None,
        "map": [
            f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
        ],
        "map_bai": [
            f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
        ],
    }
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    genome_data: NamedTuple[str | None] = lookup(
        query=f"species == '{species}' & build == '{build}' & release == '{release}'",
        within=genomes
    )

    mutect2_call_input["fasta"] = getattr(
        genome_data,
        "dna_fasta",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    )
    mutect2_call_input["fasta_fai"] = getattr(
        genome_data,
        "dna_fai",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    )
    mutect2_call_input["fasta_dict"] = getattr(
        genome_data,
        "dna_dict",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.dict",
    )

    intervals: str | None =  getattr(genome_data, "capture_kit_bed", None)
    if intervals:
        mutect2_call_input["intervals"] = intervals

    pon: str | None = getattr(genome_data, "PoN", None)
    if pon:
        mutect2_call_input["pon"] = pon

    af_only: str | None = getattr(genome_data, "germline", None)
    af_only_tbi: str | None = getattr(genome_data, "germline_tbi", None)
    if af_only and af_only_tbi:
        mutect2_call_input["germline"] = af_only
        mutect2_call_input["germline_tbi"] = af_only_tbi

    sample_data: NamedTuple[str | None] = lookup(
        query=f"species == '{species}' & build == '{build}' & release == '{release}' & sample_id == '{sample}'",
        within=samples
    )
    normal: str | None = getattr(sample_data, "corresponding_normal", None)
    if normal:
        mutect2_call_input["map"].append(
            f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{normal}.bam"
        )
        mutect2_call_input["map_bai"].append(
            f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{normal}.bam.bai"
        )

    return mutect2_call_input

def get_fair_gatk_mutect_somatic_gatk_get_pileup_summaries_input(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
) -> dict[str, str]:
    """
    Return best input files list, according to
    GATK-GetPileupSummaries snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    gatk_get_pileup_summaries_input: dict[str, str] = {
        "map":  f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "map_bai": f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
    }
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    genome_data: NamedTuple[str | None] = lookup(
        query=f"species == '{species}' & build == '{build}' & release == '{release}'",
        within=genomes
    )

    intervals: str | None =  getattr(genome_data, "capture_kit_bed", None)
    if intervals:
        gatk_get_pileup_summaries_input["intervals"] = intervals


    af_only: str | None = getattr(genome_data, "germline", None)
    af_only_tbi: str | None = getattr(genome_data, "germline_tbi", None)
    if af_only and af_only_tbi:
        gatk_get_pileup_summaries_input["variants"] = af_only
        gatk_get_pileup_summaries_input["variants_tbi"] = af_only_tbi

    return gatk_get_pileup_summaries_input


def get_fair_gatk_mutect_somatic_gatk_filter_mutect_calls_inputt(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    GATK-FilterMutectCall snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    filter_mutect_calls_input: dict[str, str] = {
        "ref": None,
        "ref_fai": None,
        "ref_dict": None,
        "aln": f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "aln_bai": f"tmp/fair_gatk_mutect_somatic/picard_replace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
        "vcf": f"tmp/fair_gatk_mutect_somatic/gatk_mutect2/{species}.{build}.{release}.{datatype}/{sample}.vcf",
        "f1r2": f"tmp/fair_gatk_mutect_somatic/gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.tar.gz",
    }
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    genome_data: NamedTuple[str | None] = lookup(
        query=f"species == '{species}' & build == '{build}' & release == '{release}'",
        within=genomes
    )

    filter_mutect_calls_input["ref"] = getattr(
        genome_data,
        "dna_fasta",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    )
    filter_mutect_calls_input["ref_fai"] = getattr(
        genome_data,
        "dna_fai",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    )
    filter_mutect_calls_input["ref_dict"] = getattr(
        genome_data,
        "dna_dict",
        f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.dict",
    )
    intervals: str | None = getattr(genome_data, "capture_kit_bed", None)
    if intervals:
        filter_mutect_calls_input["intervals"] = intervals

    af_only: str | None = genome_data.get("af_only")
    af_only_tbi: str | None = genome_data.get("af_only_tbi")
    if af_only and af_only_tbi:
        filter_mutect_calls_input["contamination"] =  f"tmp/fair_gatk_mutect_somatic/gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table"

    return filter_mutect_calls_input

def get_fair_gatk_mutect_somatic_multiqc_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    Mutect2 snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    datatype: str = str(wildcards.datatype)
    results: dict[str, list[str]] = {
        "fastp_pair_ended": collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(query="downstream_file == downstream_file", within=samples),
        ),
        "fastp_single_ended": collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
            sample=lookup(query="downstream_file != downstream_file", within=samples),
        ),
        "fastqc_pair_ended": collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query="downstream_file == downstream_file",
                within=samples,
            ),
            stream=stream_list,
        ),
        "fastqc_single_ended": collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(query="downstream_file != downstream_file", within=samples),
        ),
        "bowtie2": [],
        "samtools": [],
        "picard_qc": [],
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["bowtie2"].append(
            f"logs/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.log"
        )
        results["samtools"].append(
            f"tmp/fair_bowtie2_mapping/samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.txt"
        )
        results["picard_qc"] += multiext(
            f"tmp/fair_bowtie2_mapping/picard_create_multiple_metrics/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )
        

    return results