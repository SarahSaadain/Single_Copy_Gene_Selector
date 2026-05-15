####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def map_reads_to_scg_library_input_reads(wildcards):
    # Function to get the read file path for a given sample from the config

    # Extract species and sample from wildcards
    # They are defined in the rule as {species} and {sample}
    species = wildcards.species
    sample = wildcards.sample

    return get_path_of_sample(species, sample)

def trigger_map_reads_to_scg_library_input(wildcards):
    species = wildcards.species
    samples = get_samples_of_species(species)

    inputs = []

    for sample in samples:
        bam = f"results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.dedupped.bam"
        bai = f"{bam}.bai"
        inputs.append(bam)
        inputs.append(bai)

    return inputs

####################################################
# Snakemake rules
####################################################

_mapper = config.get("mapping", {}).get("mapper", "minimap2")

rule trigger_map_reads_to_scg_library:
    input:
        trigger_map_reads_to_scg_library_input
    output:
        temp("results/{species}/trigger_map_reads_to_scg_library.trigger"),
    message: "Mapping TE and SCG libraries to reference genomes for {wildcards.species} completed"
    shell:
        """
        touch {output}
        """

if _mapper == "minimap2":

    rule index_library_for_mapping_minimap2:
        input:
            target="results/{species}/scg_library/{species}_scg_library.fasta"
        output:
            temp("results/{species}/scg_library/{species}_scg_library.fasta.mmi"),
        log:
            "results/{species}/scg_library/{species}_scg_library_minimap2_index.log"
        message: "Indexing SCG and TE library {input} with minimap2"
        wrapper:
            "v9.3.0/bio/minimap2/index"

    rule map_reads_to_scg_library:
        input:
            query=map_reads_to_scg_library_input_reads,
            target="results/{species}/scg_library/{species}_scg_library.fasta.mmi",
        output:
            temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.with_unmapped.bam"),
        log:
            "results/{species}/reads/mapped_scg_library/{sample}_minimap2.log",
        message: "Mapping reads of {wildcards.sample} to {wildcards.species} SCG and TE library with minimap2"
        params:
            extra="-ax sr",
            sorting="coordinate",
        threads: 10
        wrapper:
            "v9.3.0/bio/minimap2/aligner"

else:
    # bwa-mem2 (fallback)

    rule index_library_for_mapping_bwa_mem2:
        input:
            "results/{species}/scg_library/{species}_scg_library.fasta"
        output:
            temp("results/{species}/scg_library/{species}_scg_library.fasta.0123"),
            temp("results/{species}/scg_library/{species}_scg_library.fasta.amb"),
            temp("results/{species}/scg_library/{species}_scg_library.fasta.ann"),
            temp("results/{species}/scg_library/{species}_scg_library.fasta.bwt.2bit.64"),
            temp("results/{species}/scg_library/{species}_scg_library.fasta.pac"),
        log:
            "results/{species}/scg_library/{species}_scg_library_bwa_index.log"
        message: "Indexing SCG and TE library {input} with BWA-MEM2"
        wrapper:
            "v9.3.0/bio/bwa-mem2/index"

    rule map_reads_to_scg_library:
        input:
            reads=map_reads_to_scg_library_input_reads,
            idx=multiext("results/{species}/scg_library/{species}_scg_library.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        output:
            temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.with_unmapped.bam"),
        log:
            "results/{species}/reads/mapped_scg_library/{sample}_bwa.log",
        message: "Mapping reads of {wildcards.sample} to {wildcards.species} SCG and TE library with BWA-MEM2"
        params:
            sort="samtools",
            sort_order="coordinate",
        threads: 10
        wrapper:
            "v9.3.0/bio/bwa-mem2/mem"


rule remove_unmapped_reads_from_bam_reads_to_library:
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.with_unmapped.bam"
    output:
        bam="results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam"
    message: "Converting SAM to BAM for {input}"
    params:
        extra="-b -F 4",  # optional params string
    threads: 2
    wrapper:
        "v9.3.0/bio/samtools/view"

# Rule: Index BAM file
# SAMTOOLS doesn't parallelize the indexing work — it only parallelizes compression/decompression.
rule index_bam_reads_to_library:
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam"
    output:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam.bai"
    message: "Indexing BAM file for {input}"
    params:
        extra="",
    threads: 5
    wrapper:
        "v9.3.0/bio/samtools/index"

