# Single_Copy_Gene_Selector
This repository contains a Snakemake pipeline which can be used to determine the best SCGs for aDNA downstream analysis.

## SCG Selection with Mixed Modern and Ancient Reads

### The Core Problem

Single-copy genes (SCGs) are widely used in ancient DNA (aDNA) analysis as calibrators for things like depth normalization, sex determination, or copy number estimation. The naïve approach — just taking whatever BUSCO calls "Complete" — turns out to be insufficient when you're working with ancient samples. Here's why, and what this pipeline does instead.

### Step 1: Identify Candidate SCGs via BUSCO (from Modern Reference)

The pipeline starts with BUSCO run against a **modern reference genome**. BUSCO identifies genes that are "Complete" and single-copy within the lineage's ortholog database (e.g., `drosophilidae_odb12`). From these, the pipeline extracts the actual nucleotide sequences of each SCG using the genomic coordinates in BUSCO's `full_table.tsv`, filtering out any that are too short (default: ≥1,500–2,000 bp). These sequences form the **SCG library** — essentially a custom miniature reference of candidate single-copy loci.

This step is necessary but not sufficient. BUSCO tells you that a gene appears single-copy in a high-quality modern genome assembly, but it says nothing about how well that gene will actually behave when you map real, degraded, ancient reads onto it.

### Step 2: Map All Samples (Modern + Ancient) to the SCG Library

Every sample — both modern reference-quality samples and ancient degraded samples — is independently mapped to this SCG library. For each sample, the pipeline computes the following statistics per SCG using `pysam`'s `count_coverage()`:

- **Min depth** — the minimum read depth at any position in the gene
- **Average depth** — mean read depth across all positions
- **Median depth** — median read depth across all positions
- **Max depth** — the peak pile-up at any single position
- **Covered bases** — the absolute number of positions with depth > 0
- **Breadth of coverage** — the proportion of the gene's length covered by at least one read (`covered_bases / length`)

Note that `count_coverage()` counts per-base observations (A/C/G/T) rather than spanning reads, so reads with deletions at a given position contribute 0 depth there. This differs subtly from `samtools depth`, but the difference is negligible for typical coverage estimation.

These stats are stored per-sample as JSON files, giving a rich, empirically-grounded picture of how each SCG actually behaves across the diversity of your dataset.

### Step 3: Score and Rank SCGs Across All Samples

This is where the pipeline diverges most sharply from a naive BUSCO-only approach. The `determine_scg_ranking.py` script aggregates stats across all BAM files and scores each SCG on three jointly penalized criteria.

**Breadth of coverage** (`score_breadth`) is the mean breadth across all samples. A gene that isn't reliably covered across the length of its sequence is useless for normalization, even if it's deeply covered in patches.

**Depth variation** (`score_depth_variation`) is penalized exponentially using the ratio of `max_depth / mean_median_depth` (the maximum depth across all samples divided by the mean of per-sample median depths). This is computed as:

```
score_depth_variation = exp(-max_variation × 0.3)
```

A true single-copy gene should have relatively uniform depth. Genes with extreme local pile-ups — caused by repetitive elements (such as microsatellites embedded within or near the gene), alignment artifacts, or paralogs missed by BUSCO — are downranked harshly.

**Depth consistency** (`score_depth_consistency`) penalizes SCGs whose depth deviates from the global median. A MAD-based (Median Absolute Deviation) deviation score is computed across all SCGs:

```
depth_deviation = |median_depth_scg - global_median_depth| / (global_MAD + ε)
score_depth_consistency = exp(-depth_deviation / 3)
```

where `global_median_depth` and `global_MAD` are derived from the distribution of per-SCG median depths (each being the mean of per-sample median depths). SCGs that are consistently under- or over-represented relative to the bulk of SCGs are likely not truly single-copy in practice, even if BUSCO said they were.

The final score is:

```
score = score_breadth + score_depth_variation + score_depth_consistency
```

SCGs are then ranked by this score in descending order and written to both a TSV and a detailed JSON summary.

### Why This Is Better Than Either Alternative Alone

**vs. BUSCO only:** BUSCO operates on a polished modern assembly under ideal conditions. It cannot know how a gene will behave under ancient DNA damage (cytosine deamination, fragmentation, short reads), how its mappability holds up with aDNA read lengths of 30–80 bp, or whether it harbours low-complexity regions that create alignment pile-ups. A gene "Complete" in BUSCO may perform terribly as a coverage calibrator in practice — for example, genes harbouring microsatellite tracts can produce severe localised pile-ups that inflate depth metrics and distort any normalization that relies on them.

Example of microsatelite in a busco gene from Clup:

![Example of microsatelite plotted with teplotter](docs/img/microsatelite_example_clup_21182at33554_NC_049239.1_46153646_46163159_scg_teplotter.png)
*Example of microsatelite plotted with teplotter*


![Example of microsatelite in IGV](docs/img/microsatelite_example_clup_21182at33554_NC_049239.1_46153646_46163159_scg_icv.png)
*Example of microsatelite in IGV*


**vs. Modern DNA only:** Modern reads are long, high-quality, and undamaged. A gene that looks clean with modern Illumina reads may behave very differently when you're aligning 50 bp, damage-saturated ancient reads. Paralogous regions that modern reads span uniquely (with longer reads or paired-end information) may cause multi-mapping disasters with short ancient reads. Using only modern samples to select SCGs creates a systematic bias where the selected genes are optimized for the wrong data type.

**The mixed approach** forces every SCG to earn its place by performing consistently across both data types simultaneously. A gene that scores well across modern and ancient samples is one that is genuinely single-copy, mappable with short reads, free of alignment artifacts under real aDNA conditions, and reliably covered at the breadths needed for downstream analysis. This makes the final SCG set much more robust as a reference for any depth-based inference in the aDNA study.