# prot2gene 🍗 ➡️ 🧬

**prot2gene** is a computational pipeline to map protein sequences to genomic coordinates, providing a simple and generalizable framework for integrative proteogenomic analysis.

The prot2gene.ipynb notebook includes code for: 

1. Extracting protein-coding regions from the human reference genome annotation file (GTF)
2. Generating a dataset of protein translations directly from coding sequences (CDS) in the reference genome
3. Tracking protein sequences to genomic CDS coordinates

**Example usage:**

```python
query_protein = "DLVILLYETALLSSGFSLEDPQTHANRIYRMIKLGLGIDEDDPTADDTSAAVTEEMPPLE" # input any protein sequence here

gencode_cds_coordinates["matched_cds_coords"] = gencode_cds_coordinates.apply(
    lambda row: find_cds_coords_for_protein_subseq(row, query_protein),
    axis=1
)

matches = gencode_cds_coordinates[gencode_cds_coordinates["matched_cds_coords"].notnull()]

matches[["gene_name", "seqname", "strand", "gene_id", "transcript_id", "protein_id", "matched_cds_coords"]]
```
**Output:**
```
	gene_name	seqname	strand	gene_id	transcript_id	protein_id	matched_cds_coords
1387198	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000216281.13	ENSP00000216281.8	[[[102081751, 102081821], [102082111, 10208221...
1387224	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000334701.11	ENSP00000335153.7	[[[102081751, 102081821], [102082111, 10208221...
```

Below is a detailed explanation of the prot2gene pipeline implementation.

# Extracting protein-coding regions from the human reference genome annotation file (GTF)

