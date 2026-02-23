# Single_Copy_Gene_Selector
This repository contains a Snakemake pipeline which can be used to determine the best SCGs for aDNA downstream Analysis

Now I have a thorough understanding of the pipeline. Let me explain the process clearly.

## SCG Selection with Mixed Modern and Ancient Reads

### The Core Problem

Single-copy genes (SCGs) are widely used in ancient DNA (aDNA) analysis as calibrators for things like depth normalization, sex determination, or copy number estimation. The naïve approach — just taking whatever BUSCO calls "Complete" — turns out to be insufficient when you're working with ancient samples. Here's why, and what this pipeline does instead.

---

### Step 1: Identify Candidate SCGs via BUSCO (from Modern Reference)

The pipeline starts with BUSCO run against a **modern reference genome**. BUSCO identifies genes that are "Complete" and single-copy within the lineage's ortholog database (e.g., `drosophilidae_odb12`). From these, the pipeline extracts the actual nucleotide sequences of each SCG using the genomic coordinates in BUSCO's `full_table.tsv`, filtering out any that are too short (default: ≥1,500–2,000 bp). These sequences form the **SCG library** — essentially a custom miniature reference of candidate single-copy loci.

This step is necessary but not sufficient. BUSCO tells you that a gene appears single-copy in a high-quality modern genome assembly, but it says nothing about how well that gene will actually behave when you map real, degraded, ancient reads onto it.

---

### Step 2: Map All Samples (Modern + Ancient) to the SCG Library

Every sample — both modern reference-quality samples and ancient degraded samples — is independently mapped to this SCG library. For each sample, the pipeline computes three statistics per SCG:

- **Average depth** — how many reads cover the gene on average
- **Maximum depth** — the peak pile-up at any single position
- **Breadth of coverage** — the proportion of the gene's length covered by at least one read

These are stored per-sample as JSON files, giving a rich, empirically-grounded picture of how each SCG actually behaves across the diversity of your dataset.

---

### Step 3: Score and Rank SCGs Across All Samples

This is where the pipeline diverges most sharply from a naive BUSCO-only approach. The `determine_scg_ranking.py` script aggregates stats across all BAM files and scores each SCG on three jointly penalized criteria:

**Breadth of coverage** is weighted heavily — a gene that isn't reliably covered across the length of its sequence is useless for normalization, even if it's deeply covered in patches.

**Depth variation** is penalized exponentially. A true single-copy gene should have relatively uniform depth. Genes with extreme local pile-ups (often caused by repetitive elements, alignment artifacts, or paralogs missed by BUSCO) are downranked harshly. The ratio of `max_depth / mean_depth` catches these problematic genes.

**Depth deviation from the median SCG** is penalized using a MAD-based (Median Absolute Deviation) z-score — a robust, outlier-resistant statistic. SCGs that are consistently under- or over-represented relative to the bulk of SCGs are likely not truly single-copy in practice, even if BUSCO said they were. This catches cases like recently duplicated genes or genes in difficult genomic regions.

The final score is: `breadth + exp(-variation) + exp(-depth_deviation)`. The top-N SCGs (e.g., 100) are then filtered to a final ranked list.

---

### Why This Is Better Than Either Alternative Alone

**vs. BUSCO only:** BUSCO operates on a polished modern assembly under ideal conditions. It cannot know how a gene will behave under ancient DNA damage (cytosine deamination, fragmentation, short reads), how its mappability holds up with aDNA read lengths of 30–80 bp, or whether it harbors low-complexity regions that create alignment pile-ups. A gene "Complete" in BUSCO may perform terribly as a coverage calibrator in practice.

**vs. Modern DNA only:** Modern reads are long, high-quality, and undamaged. A gene that looks clean with modern Illumina reads may behave very differently when you're aligning 50 bp, damage-saturated ancient reads. Paralogous regions that modern reads span uniquely (with longer reads or paired-end information) may cause multi-mapping disasters with short ancient reads. Using only modern samples to select SCGs creates a systematic bias where the selected genes are optimized for the wrong data type.

**The mixed approach** forces every SCG to earn its place by performing consistently across both data types simultaneously. A gene that scores well across modern and ancient samples is one that is genuinely single-copy, mappable with short reads, free of alignment artifacts under real aDNA conditions, and reliably covered at the breadths needed for downstream analysis. This makes the final SCG set much more robust as a reference for any depth-based inference in the aDNA study.