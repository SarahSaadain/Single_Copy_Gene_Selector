import pysam
import collections
import logging
import sys
import json

bam_file = snakemake.input.bam  # single BAM
log_filename = snakemake.log[0]
stats_file = snakemake.output.stats  # intermediate JSON with SCG stats

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stderr), logging.FileHandler(log_filename)]
)
logger = logging.getLogger(__name__)

def get_coverage_stats(bam_path, contig, length):
    """Compute mean depth, max depth, and breadth for a contig in a BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        depth_sum = 0
        depth_max = 0
        covered_positions = 0
        for pileupcolumn in bam.pileup(contig, start=0, stop=length, truncate=True, stepper="all"):
            depth = pileupcolumn.nsegments
            depth_sum += depth
            if depth > 0:
                covered_positions += 1
            if depth > depth_max:
                depth_max = depth
        avg_depth = depth_sum / length if length > 0 else 0.0
        breadth = covered_positions / length if length > 0 else 0.0
        return avg_depth, depth_max, breadth

# collect SCG stats
scg_stats = {}

logger.info(f"Processing BAM: {bam_file}")
with pysam.AlignmentFile(bam_file, "rb") as bam:
    for name, length in zip(bam.references, bam.lengths):
        avg_depth, max_depth, breadth = get_coverage_stats(bam_file, name, length)
        scg_stats[name] = {
            "avg_depth": avg_depth,
            "max_depth": max_depth,
            "breadth": breadth
        }

# save to JSON
with open(stats_file, "w") as f:
    json.dump(scg_stats, f, indent=2)

logger.info(f"SCG stats for {bam_file} saved to {stats_file}")
