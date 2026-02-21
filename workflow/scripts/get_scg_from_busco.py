import os
from typing import Tuple, Optional, List
import logging
from Bio.Seq import Seq
import pysam
import pandas as pd
from snakemake.script import snakemake


def extract_busco_sequences(df, reference_fasta, min_length_scg) -> List[Tuple[str, str]]:
    fasta_entries = []

    processed_seqs = 0

    for i, row in df.iterrows():

        try:
            busco_id = row["Busco_id"]
            chrom = row["Sequence"]
            start = int(row["Gene_Start"])
            end = int(row["Gene_End"])
            strand = row["Strand"]

            # Skip rows missing key fields
            if pd.isna(chrom) or pd.isna(start) or pd.isna(end):
                logger.warning(f"Skipping {busco_id}: missing coordinates.")
                continue

            # if compliment strand, swap start and end
            if strand == "-":
                start, end = end, start

            if end - start <= min_length_scg:
                logger.warning(f"Skipping {busco_id}: length {end - start} is less than minimum {min_length_scg}.")
                continue

            seq = extract_sequence(reference_fasta, chrom, start, end)
            header = f"{busco_id}_{chrom}_{start}_{end}_scg"
            fasta_entries.append((header, seq))

            processed_seqs += 1

        except Exception as e:
            logger.error(f"Error processing {busco_id}: {e}")
    
    return fasta_entries

def extract_sequence(reference_fasta: str, chrom: str, start: int, end: int) -> str:
    fasta = pysam.FastaFile(reference_fasta)
    try:
        sequence = fasta.fetch(chrom, start, end)
    except Exception as e:
        raise RuntimeError(f"Could not extract {chrom}:{start}-{end} from FASTA: {e}")
    
    return sequence


def write_multi_fasta(out_file: str, entries: List[Tuple[str, str]]):

    with open(out_file, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n{seq}\n")


def main():

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stderr)]
    )
    logger = logging.getLogger(__name__)

    full_table_path = snakemake.input.busco_full_table
    reference_fasta = snakemake.input.ref_genome[0]
    min_length_scg = snakemake.params.min_length_scg
    out_file = snakemake.output.scg

    if not os.path.isfile(full_table_path):
        raise ValueError(f"full_table.tsv not found in {full_table_path}")

    if not os.path.isfile(reference_fasta):
        raise ValueError(f"Reference FASTA not found: {reference_fasta}")
    
    
    # load full_table.tsv into a pandas dataframe and skip first 3 lines
    # there are no headers
    # Sample:
    # # BUSCO version is: 6.0.0 
    # # The lineage dataset is: drosophilidae_odb12 (Creation date: 2025-07-01, number of genomes: 8, number of BUSCOs: 9264)
    # # Busco id      Status  Sequence        Gene Start      Gene End        Strand  Score   Length  OrthoDB url     Description
    # 11at7214        Complete        JAEIFK010000455.1       869678  869976  +       211.1   94      https://www.orthodb.org/?query=11at7214 cytochrome c oxidase subunit 6B1
    # 28at7214        Complete        JAEIFK010000678.1       13548246        13548860        +       151.2   55      https://www.orthodb.org/?query=28at7214 Small ribosomal subunit protein uS14
    # 31at7214        Complete        JAEIFK010000425.1       9934295 9935384 +       196.3   81      https://www.orthodb.org/?query=31at7214 protein drumstick
    # 97at7214        Complete        JAEIFK010000675.1       1798636 1805396 +       968.7   435     https://www.orthodb.org/?query=97at7214 neurogenic locus notch homolog protein 1

        
    # read table while skipping comment lines
    df = pd.read_csv(
        full_table_path,
        sep='\t',
        comment='#',
        header=None,
        names=[
            "Busco_id", "Status", "Sequence", "Gene_Start", "Gene_End",
            "Strand", "Score", "Length", "OrthoDB_url", "Description"
        ],
        dtype=str  # keep as string to avoid mixed types
    )

    # Keep only complete BUSCOs with valid coordinates
    df = df[df["Status"].str.contains("Complete", na=False)]

    fasta_entries = extract_busco_sequences(df, reference_fasta, min_length_scg)

    logger.info(f"Extracted {len(fasta_entries)} SCG sequences from BUSCO results.")

    write_multi_fasta(out_file, fasta_entries)
    logger.info(f"Sequences written to: {out_file}")


if __name__ == "__main__":
    main()
