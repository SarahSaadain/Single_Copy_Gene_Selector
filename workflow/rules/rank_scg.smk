rule determine_scg_ranking_across_samples:
    input:
        stats = lambda wildcards: expand("results/{species}/scg_library/{sample}_scg_stats.json",
                species=wildcards.species,
                sample=get_samples_of_species(wildcards.species)
            )
    output:
        best_scgs = "results/{species}/{species}_scg_ranked.tsv",
        best_scgs_json = "results/{species}/{species}_scg_ranked.json"
    message: "Determining best SCG for {input}"
    log: "results/{species}/{species}_scg_ranked.log"
    conda:
        "../envs/python.yaml"
    script: "../scripts/determine_scg_ranking.py"

rule compute_scg_stats_for_bam:
    input:
        bam = "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam",
        bai = "results/{species}/reads/mapped_scg_library/{sample}_scg_library.sorted.bam.bai"
    output:
        stats = "results/{species}/scg_library/{sample}_scg_stats.json"
    message: "Determining best SCG for {input}"
    log: "results/{species}/scg_library/{sample}_scg_stats.log"
    conda:
        "../envs/python.yaml"
    script: "../scripts/compute_scg_stats_for_bam.py"

rule filter_top_scgs_tes:
    input:
        ranked_scgs="results/{species}/{species}_scg_ranked.tsv"
    output:
        temp_picked_scgs=temp("results/{species}/{species}_top_scgs.tsv"),
        relevant_contigs="results/{species}/{species}_relevant_scg.txt",
        relevant_contigs_bed="results/{species}/{species}_relevant_scg.bed"
    params:
        num_top_scgs=lambda wildcards: config["species"][wildcards.species].get("settings", {}).get("num_top_scgs", 50)
    shell:
        """
        awk 'NR>1 && NR<={params.num_top_scgs}+1 {{print $1}}' {input.ranked_scgs} | sort -u > {output.relevant_contigs}
        cp {output.relevant_contigs} {output.temp_picked_scgs}
        awk '{{print $1 "\t0\t1000000000"}}' {output.relevant_contigs} | sort -u > {output.relevant_contigs_bed}
        """

rule filter_fasta:
    input:
        fasta = "results/{species}/scg_library/{species}_scg_library.fasta",
        id_list = "results/{species}/{species}_relevant_scg.txt"
    output:
        filtered = "results/{species}/{species}_relevant_scg.fasta"
    run:
        # Read the list of IDs to keep
        with open(input.id_list) as f:
            ids_to_keep = set(line.strip() for line in f if line.strip())

        # Parse and filter the FASTA
        with open(input.fasta) as fin, open(output.filtered, "w") as fout:
            write = False
            for line in fin:
                if line.startswith(">"):
                    # Match just the ID (before first space)
                    seq_id = line[1:].split()[0].strip()
                    write = seq_id in ids_to_keep
                if write:
                    fout.write(line)
