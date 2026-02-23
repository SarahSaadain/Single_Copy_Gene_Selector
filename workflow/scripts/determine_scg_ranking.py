import json
import numpy as np
import logging
import sys

stats_files = snakemake.input.stats  # list of JSON files from Step 1
summary_file = snakemake.output.best_scgs
log_filename = snakemake.log[0]

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stderr), logging.FileHandler(log_filename)]
)
logger = logging.getLogger(__name__)

# Load all SCG stats
scg_stats_per_bam = []  # list of dicts
for stats_file in stats_files:
    with open(stats_file) as f:
        scg_stats_per_bam.append(json.load(f))

# Aggregate across BAMs
from collections import defaultdict

scg_aggregated = defaultdict(list)

for scg_stats in scg_stats_per_bam:
    for scg, stats in scg_stats.items():
        scg_aggregated[scg].append(stats)

# Compute mean metrics across BAMs
scg_summary = {}
for scg, values in scg_aggregated.items():
    avg_depths = [v["avg_depth"] for v in values]
    max_depths = [v["max_depth"] for v in values]
    breadths = [v["breadth"] for v in values]
    mean_depth = np.mean(avg_depths)
    max_variation = max(max_depths) / mean_depth if mean_depth > 0 else 0
    mean_breadth = np.mean(breadths)
    scg_summary[scg] = {
        "mean_depth": mean_depth,
        "max_variation": max_variation,
        "mean_breadth": mean_breadth
    }

# Normalize metrics for scoring
breadths = [v["mean_breadth"] for v in scg_summary.values()]
variations = [v["max_variation"] for v in scg_summary.values()]
min_var, max_var = min(variations), max(variations)

all_mean_depths = [v["mean_depth"] for v in scg_summary.values()]
median_depth = np.median(all_mean_depths)
mad_depth = np.median(np.abs(np.array(all_mean_depths) - median_depth))  # median absolute deviation, robust to outliers

scg_scores = {}
for scg, vals in scg_summary.items():
    norm_breadth = vals["mean_breadth"]  # weight breadth more heavily, tune as needed
    
    # Exponential penalty for depth variation
    mean_var = np.mean(variations)
    norm_var = np.exp(-vals["max_variation"] / mean_var) if mean_var > 0 else 1
    
    # Penalty for depth being far from the median SCG depth
    # Using MAD-based z-score (robust to outliers)
    depth_deviation = abs(vals["mean_depth"] - median_depth) / (mad_depth + 1e-9)
    norm_depth = np.exp(-depth_deviation / 2)  # tune the /2 to be stricter or looser
    
    score = norm_breadth + norm_var + norm_depth
    scg_scores[scg] = score

# Select top N SCGs
top_scgs = sorted(scg_scores.items(), key=lambda x: x[1], reverse=True)

# Write summary
with open(summary_file, "w") as f:
    f.write("scg_name\tmean_depth\tmax_depth_variation\tmean_breadth\tscore\n")
    for scg, score in top_scgs:
        stats = scg_summary[scg]
        f.write(f"{scg}\t{stats['mean_depth']:.2f}\t{stats['max_variation']:.2f}\t{stats['mean_breadth']:.3f}\t{score:.3f}\n")
