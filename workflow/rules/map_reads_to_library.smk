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


# This rule indexes the combined SCG and TE library using BWA2 for read mapping
rule index_library_for_mapping:
    input:
        "results/{species}/scg_library/{species}_scg_library.fasta"
    output:
        "results/{species}/scg_library/{species}_scg_library.fasta.0123",
        "results/{species}/scg_library/{species}_scg_library.fasta.amb",
        "results/{species}/scg_library/{species}_scg_library.fasta.ann",
        "results/{species}/scg_library/{species}_scg_library.fasta.bwt.2bit.64",
        "results/{species}/scg_library/{species}_scg_library.fasta.pac",
    log:
        "results/{species}/scg_library/{species}_scg_library_bwa_index.log"
    message: "Indexing SCG and TE library {input} with BWA2"
    wrapper:
        "v9.3.0/bio/bwa-mem2/index"

# This rule maps sequencing reads of each sample to the combined SCG and TE library
rule map_reads_to_scg_library:
    input:
        # Since we do not know where the reads are physically located, 
        # we use a function to get the correct path from the config file
        reads=map_reads_to_scg_library_input_reads,
        # The index is produced by the bwa_index_library rule
        idx=multiext("results/{species}/scg_library/{species}_scg_library.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    output:
        temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.sam"), # became temp to save space
    log:
        "results/{species}/reads/mapped_scg_library/{sample}_bwa.log",
    message: "Mapping reads of {wildcards.sample} to {wildcards.species} SCG and TE library"
    params:
        #extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    threads: 10
    wrapper:
        "v9.3.0/bio/bwa-mem2/mem"


rule convert_sam_to_bam_reads_to_library:
    # 2 Convert SAM to BAM
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sam",
    output:
        bam=temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.bam") # needs to be called .unsorted.bam otherwise snakemake had problems with unambigous names
    message: "Converting SAM to BAM for {input}"
    threads: 8
    wrapper:
        "v9.3.0/bio/samtools/view"

rule remove_unmapped_reads_from_bam_reads_to_library:
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.bam"
    output:
        bam=temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.no_unmapped.bam") # needs to be called .unsorted.bam otherwise snakemake had problems with unambigous names
    message: "Converting SAM to BAM for {input}"
    params:
        extra="-b -F 4",  # optional params string
    threads: 2
    wrapper:
        "v9.3.0/bio/samtools/view"

rule  sort_bam_reads_to_library:
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.no_unmapped.bam"
    output:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam"
    message: "Sorting BAM file for {input}"
    log:
        "results/{species}/reads/mapped_scg_library/{sample}_sort_bam.log",
    threads: 8
    wrapper:
        "v9.3.0/bio/samtools/sort"

rule deduplicate_bam_with_dedup:
    input:
        bam         = "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam"
    output:
        dedup_folder= directory("results/{species}/reads/mapped_scg_library/deduplication/{sample}/"),
        dedup_bam   = "results/{species}/reads/mapped_scg_library/deduplication/{sample}/{sample}_scg_library.sorted_rmdup.bam",
        dedup_hist  = "results/{species}/reads/mapped_scg_library/deduplication/{sample}/{sample}_scg_library.sorted.hist",
        dedup_json  = "results/{species}/reads/mapped_scg_library/deduplication/{sample}/{sample}_scg_library.sorted.dedup.json",
    message:
        "Deduplicating BAM file for {input.bam} using dedup for sample {wildcards.sample} in species {wildcards.species}",
    log: 
        "results/{species}/reads/mapped_scg_library/deduplication/{sample}_dedup.log",
    conda:
        "../envs/dedup.yaml"
    resources:
        mem_mb = 20000   # request 10 GB from cluster / cgroups
    shell:
        """
        mkdir -p {output.dedup_folder}
        # Set explicit heap size via -Xms (initial) and -Xmx (max)
        dedup -Xms5g -Xmx20g --input {input.bam} --merged --output {output.dedup_folder} > {log}
        """

rule move_deduplicated_to_library_mapped:
    input:
        dedup_bam = "results/{species}/reads/mapped_scg_library/deduplication/{sample}/{sample}_scg_library.sorted_rmdup.bam",
    output:
        moved_bam = temp("results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.dedupped.bam")
    shell:
        """
        mv {input.dedup_bam} {output.moved_bam}
        """

# Rule: Sort BAM file
rule sort_dedup_bam_reads_to_library:
    # 3 Sort BAM
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.unsorted.dedupped.bam"
    output:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.dedupped.bam"
    message: "Sorting BAM file for {input}"
    log:
        "results/{species}/reads/mapped_scg_library/{sample}.sorted.dedupped.bam.log",
    threads: 10
    wrapper:
        "v9.3.0/bio/samtools/sort"

# Rule: Index BAM file
# SAMTOOLS doesn’t parallelize the indexing work — it only parallelizes compression/decompression.
rule index_bam_reads_to_library:
    # 4 Index BAM
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam" 
    output:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam.bai"
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v9.3.0/bio/samtools/index"

# Rule: Index BAM file
# SAMTOOLS doesn’t parallelize the indexing work — it only parallelizes compression/decompression.
rule index_dedup_bam_reads_to_library:
    # 4 Index BAM
    input:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.dedupped.bam" 
    output:
        "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.dedupped.bam.bai"
    message: "Indexing BAM file for {input}"
    params: 
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v9.3.0/bio/samtools/index"
