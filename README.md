# prot2gene 🍗 ➡️ 🧬

**prot2gene** is a computational pipeline to map protein sequences to genomic coordinates, providing a simple and generalizable framework for integrative proteogenomic analysis.

This repository includes a **pre-computed CDS-to-protein database for all protein-coding transcripts in the human genome**: gencode_cds_coordinates.zip 

This file can be used directly to map protein sequences to genomic coordinates in **2 steps**: 

1. Define the **prot2gene function**:

```python

def prot2gene(row, query_protein):
    """
    Given a row with 'cds_sequence', 'protein_sequence', and 'cds_coordinates',
    find the genomic CDS coordinates corresponding to a peptide subsequence.
    Returns a list of coordinate lists (one per match).
    """
    protein_seq = row['protein_sequence']
    cds_seq = row['cds_sequence']
    strand = row['strand']
    coords = row['cds_coordinates']

    if not protein_seq or not cds_seq or not coords:
        return None

    # Find all occurrences of the query peptide
    matches = []
    start_idx = 0
    while True:
        idx = protein_seq.find(query_protein, start_idx)
        if idx == -1:
            break
        matches.append(idx)
        start_idx = idx + 1

    if not matches:
        return None

    # Reverse order of CDS exons if strand is negative
    exon_coords = coords[::-1] if strand == '-' else coords

    result_coords = []
    cds_len = sum(end - start + 1 for start, end in exon_coords)

    for prot_start in matches:
        prot_end = prot_start + len(query_protein)
        nuc_start = prot_start * 3
        nuc_end = prot_end * 3

        sub_coords = []
        cum_length = 0

        for start, end in exon_coords:
            exon_len = end - start + 1
            exon_start_in_cds = cum_length
            exon_end_in_cds = cum_length + exon_len

            # Determine overlap with the peptide's CDS interval
            overlap_start = max(exon_start_in_cds, nuc_start)
            overlap_end = min(exon_end_in_cds, nuc_end)

            if overlap_start < overlap_end:
                # Offset within exon
                offset_start = overlap_start - exon_start_in_cds
                offset_end = overlap_end - exon_start_in_cds

                if strand == '+':
                    genome_start = start + offset_start
                    genome_end = start + offset_end - 1
                else:
                    # For negative strand, reverse mapping
                    genome_end = end - offset_start
                    genome_start = end - offset_end + 1

                sub_coords.append([genome_start, genome_end])

            cum_length += exon_len

        # Reorder back to genomic order
        if strand == '-':
            sub_coords = sub_coords[::-1]

        result_coords.append(sub_coords)

    return result_coords

```

2. Load the **pre-computed dataset** using pandas:

```python

import pandas as pd

gencode_cds_coordinates = pd.read_csv("gencode_cds_coordinates.zip", compression="zip")

```

Once the **prot2gene function** is defined and **gencode_cds_coordinates dataframe** is loaded, **prot2gene can be used to map the genomic coordinates for any protein sequence in the human genome**. 

**Example usage:**

```python
query_protein = "DLVILLYETALLSSGFSLEDPQTHANRIYRMIKLGLGIDEDDPTADDTSAAVTEEMPPLE" # input any protein sequence here

gencode_cds_coordinates["matched_cds_coords"] = gencode_cds_coordinates.apply( 
    lambda row: prot2gene(row, query_protein), 
    axis=1
)                                                                              # apply find_cds_coords_for_protein function to gencode_cds_coordinates

matches = gencode_cds_coordinates[gencode_cds_coordinates["matched_cds_coords"].notnull()]

matches[["gene_name", "seqname", "strand", "gene_id", "transcript_id", "protein_id", "matched_cds_coords"]]
```
**Output:**

```
	gene_name	seqname	strand	gene_id	transcript_id	protein_id	matched_cds_coords
1387198	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000216281.13	ENSP00000216281.8	[[[102081751, 102081821], [102082111, 10208221...
1387224	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000334701.11	ENSP00000335153.7	[[[102081751, 102081821], [102082111, 10208221...
```

**The provided function and pre-computed dataset work for identifying the genomic coordinates for any annotated protein sequence from the human reference genome (GENCODE v44).** 

For additional datasets, the prot2gene.ipynb notebook includes code for: 

1. Extracting protein-coding regions from the human reference genome annotation file (GTF)
2. Generating a dataset of protein translations directly from coding sequences (CDS) in the reference genome
3. Tracking protein sequences to genomic CDS coordinates

This method provides a simple and generalizable framework for mapping protein sequences to genomic coordinates across multiple datasets. 

Below is step-by-step description of the workflow used to build the pre-computed dataset from the human reference genome.

# Extracting protein-coding regions from the human reference genome annotation file (GTF)

