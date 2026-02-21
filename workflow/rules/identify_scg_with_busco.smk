####################################################
# Snakemake rules
####################################################

rule prepare_reference:
    input:
        ref=lambda wildcards: config["species"][wildcards.species]["reference"]
    output:
        ref_link="results/{species}/ref/{reference}.fasta"
    message: "Preparing reference genome {wildcards.reference} for {wildcards.species}"
    shell:
        """
        ln -s {input.ref} {output.ref_link}
        """

rule run_busco:
    input:
        lambda wildcards: 
            expand(
                "results/{species}/ref/{reference}.fasta",
                species = wildcards.species,
                reference = get_reference_id_of_species(wildcards.species)
            )
    output:
        short_json="results/{species}/busco/short_summary.json",
        short_txt="results/{species}/busco/short_summary.txt",
        full_table="results/{species}/busco/full_table.tsv",
        miss_list="results/{species}/busco/busco_missing.tsv",
        out_dir = temp(directory("results/{species}/busco/output/")),
        dataset_dir=temp(directory("results/{species}/busco/busco_downloads")),
    log:
        "results/{species}/busco/busco.log",
    message: "Running BUSCO to identify single copy genes for {wildcards.species}"
    params:
        mode="genome",
        lineage=lambda wildcards: config["species"][wildcards.species]["lineage"],
        # optional parameters
        extra="",
    threads: 8
    wrapper:
        "v7.3.0/bio/busco"

rule prepare_scg_library:
    input:
        # Get the BUSCO full table for the species. This file contains information about
        # the identified single-copy genes (SCG) from the BUSCO analysis including their 
        # location in the reference genome
        busco_full_table="results/{species}/busco/full_table.tsv",
        # Get the first reference genome for the species from the config file.
        # Using the location from Busco, we can extract the sequences of the SCG from this genome
        ref_genome=lambda wildcards: 
            expand(
                "results/{species}/ref/{reference}.fasta",
                species = wildcards.species,
                reference = get_reference_id_of_species(wildcards.species)
            )
    output:
        scg="results/{species}/scg_library/{species}_scg_library.fasta"
    message: "Getting SCGs for {wildcards.species}"
    conda:
        "../envs/python.yaml"
    params:
        min_length_scg=lambda wildcards: config["species"][wildcards.species].get("settings", {}).get("min_length_scg", 2000)
    script:
        "../scripts/get_scg_from_busco.py"