# Pipeline Setup Guide

## Requirements

The only manual prerequisite is [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0. All other dependencies (BUSCO, BWA, SAMtools, pysam, numpy, etc.) are automatically installed by the pipeline when it runs.

To install Snakemake:

```bash
conda install -c bioconda -c conda-forge snakemake
```

## Input Data

The pipeline requires two inputs per species:

**Reference genome** — a FASTA file of the assembly against which BUSCO will be run to identify candidate SCGs. This should be a high-quality modern genome assembly.

**Reads** — one or more FASTQ files (gzipped) containing the reads to be mapped to the SCG library. These can be a mix of modern and ancient samples; providing both is strongly recommended (see the [main README](README.md) for why).

## Configuration

All pipeline parameters are controlled through `config.yaml`. A minimal working example is shown below:

```yaml
# config.yaml - Configuration file for SCG Selector Workflow
species:
  demo:
    name: demo
    lineage: drosophilidae_odb12
    settings:
      num_top_scgs: 100
      min_length_scg: 1500
    reads:
      - /demo_data/demo.fastq.gz
    reference: /demo_data/Dfun_reference_genome.fasta
```

### Configuration fields

| Field | Description |
|---|---|
| `species.<id>.name` | A short identifier for the species or run. Used to name output files. |
| `species.<id>.lineage` | The BUSCO lineage dataset to use (e.g. `drosophilidae_odb12`). Must match a dataset available in your BUSCO installation. |
| `settings.num_top_scgs` | Number of top-ranked SCGs to retain after scoring. Default: `100`. |
| `settings.min_length_scg` | Minimum SCG sequence length in base pairs. Shorter candidates are discarded. Default: `1500`. |
| `reads` | List of paths to input FASTQ files (gzipped). At least one file is required. Add one path per line for multiple samples. |
| `reference` | Path to the reference genome FASTA file. |

### Adding multiple species or samples

Multiple species blocks can coexist in the same config file, and multiple read files can be listed under a single species. Each read file is treated as an independent sample and mapped separately to the SCG library:

```yaml
species:
  drosophila:
    name: drosophila
    lineage: drosophilidae_odb12
    settings:
      num_top_scgs: 100
      min_length_scg: 1500
    reads:
      - /data/modern_sample.fastq.gz
      - /data/ancient_sample_1.fastq.gz
      - /data/ancient_sample_2.fastq.gz
    reference: /data/reference.fasta
```

## Running the Pipeline

Once your `config.yaml` is set up, run the pipeline from the repository root with:

```bash
snakemake --cores <N> --use-conda
```

Replace `<N>` with the number of CPU cores to use. For a dry run (to check that the workflow is correctly configured without executing anything):

```bash
snakemake --configfile config.yaml --cores 1 --use-conda --dry-run
```

## Output Files

The pipeline produces the following outputs per species:

| File | Description |
|---|---|
| `results/<name>/scg_stats/<sample>.json` | Per-contig coverage statistics for each input BAM (depth, breadth, etc.). |
| `results/<name>/<species>_best_scgs.tsv` | Tab-separated ranked list of SCGs with scores for breadth, depth variation, and depth consistency. |
| `results/<name>/<species>_scg_summary.json` | Full JSON summary including per-sample raw stats, aggregated metrics, and scoring details for every SCG. |
| `results/<name>/<species>_relevant_scg.fasta` | FASTA file containing the sequences of the top-ranked SCGs. |

The TSV and JSON outputs are both sorted by final score in descending order, so the top entries are the highest-quality SCGs recommended for downstream use.