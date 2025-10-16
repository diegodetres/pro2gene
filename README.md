# pro2gene 

**pro2gene** is a computational pipeline to map protein sequences to genomic coordinates, providing a simple and generalizable framework for integrative proteogenomic analysis.

This repository includes a **pre-computed CDS-to-protein database for all protein-coding transcripts in the human genome**: gencode_cds_coordinates.zip 

This file can be used directly to map protein sequences to genomic coordinates in **2 steps**: 

1. Define the **pro2gene function**:

```python

def pro2gene(row, query_protein):
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

Once the **pro2gene function** is defined and **gencode_cds_coordinates dataframe** is loaded, **pro2gene can be used to map the genomic coordinates for any protein sequence in the human genome**. 

**Example usage:**

```python
query_protein = "DLVILLYETALLSSGFSLEDPQTHANRIYRMIKLGLGIDEDDPTADDTSAAVTEEMPPLE" # input any protein sequence here

gencode_cds_coordinates["matched_cds_coords"] = gencode_cds_coordinates.apply( 
    lambda row: pro2gene(row, query_protein), 
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

For additional datasets, the pro2gene.ipynb notebook includes code for: 

1. Extracting protein-coding regions from the human reference genome annotation file (GTF)
2. Generating a dataset of protein translations directly from coding sequences (CDS) in the reference genome
3. Mapping protein sequences to genomic CDS coordinates

This method provides a simple and generalizable framework for mapping protein sequences to genomic coordinates across multiple datasets. 

Below is step-by-step description of the workflow used to build the pre-computed dataset from the human reference genome.

# Extracting protein-coding regions from the human reference genome annotation file (GTF)

Protein-coding transcripts were filtered from the **human reference genome annotation file (GTF)**, obtained from GENCODE (Release 44):

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz 
```

This file was used as a pandas dataframe and filtered for relevant transcripts:

``` python

gtf_file = 'gencode.v44.basic.annotation.gtf.gz' # GTF file

columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#', names=columns)

filtered_df = df[df['feature'] == 'transcript']
filtered_df['gene_id'] = filtered_df['attribute'].str.extract('gene_id "([^"]+)"')
filtered_df['transcript_id'] = filtered_df['attribute'].str.extract('transcript_id "([^"]+)"')
filtered_df['protein_id'] = filtered_df['attribute'].str.extract('protein_id "([^"]+)"')
filtered_df['gene_name'] = filtered_df['attribute'].str.extract('gene_name "([^"]+)"')
filtered_df = filtered_df[filtered_df['attribute'].str.contains('gene_type "protein_coding"')]

gencode_transcripts = filtered_df

```

# Generating a dataset of protein translations directly from coding sequences (CDS) in the reference genome

A gffutils SQLite database was created from the GTF file:

```python

database_filename = 'gencode_v44'

gffutils.create_db(gtf_file, database_filename)

```
This database was used to extract coding sequence (CDS) coordinates from human transcripts as follows: 

```python

db = gffutils.FeatureDB(database_filename)

df = gencode_transcripts.copy()

df['cds_coordinates'] = None

for idx, row in df.iterrows():
    stable_id = row['transcript_id']
    try:
        cds = list(db.children(stable_id, order_by='+end', featuretype=['CDS']))
        start_end_cds = [[i.start, i.end] for i in cds]
    except Exception:
        start_end_cds = []
    
    df.at[idx, 'cds_coordinates'] = start_end_cds

gencode_cds_coordinates = df

```
Sequences for the human reference genome sequence were obtained from GENCODE:

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fna.gz 
```
This file was used to map coding sequences to human transcripts as follows:

```python

fasta_path = "GCF_000001405.26_GRCh38_genomic.fna.gz" # fasta file

with gzip.open(fasta_path, "rt") as handle: 
    genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

# Map chromosome numbers to refseq IDs

chr_to_refseq = {
    "chr1": "NC_000001.11",
    "chr2": "NC_000002.12",
    "chr3": "NC_000003.12",
    "chr4": "NC_000004.12",
    "chr5": "NC_000005.10",
    "chr6": "NC_000006.12",
    "chr7": "NC_000007.14",
    "chr8": "NC_000008.11",
    "chr9": "NC_000009.12",
    "chr10": "NC_000010.11",
    "chr11": "NC_000011.10",
    "chr12": "NC_000012.12",
    "chr13": "NC_000013.11",
    "chr14": "NC_000014.9",
    "chr15": "NC_000015.10",
    "chr16": "NC_000016.10",
    "chr17": "NC_000017.11",
    "chr18": "NC_000018.10",
    "chr19": "NC_000019.10",
    "chr20": "NC_000020.11",
    "chr21": "NC_000021.9",
    "chr22": "NC_000022.11",
    "chrX": "NC_000023.11",
    "chrY": "NC_000024.10",
    "chrM": "NC_012920.1" 
}

gencode_cds_coordinates['refseq_id'] = gencode_cds_coordinates['seqname'].map(chr_to_refseq)

# Extract and concatenate CDS segments

def get_cds_sequence(row):
    chrom = row['refseq_id']
    strand = row['strand']
    coords = row['cds_coordinates']
    
    if chrom not in genome or not coords:
        return None
    
    # Concatenate all CDS segments
    seq = ''.join(str(genome[chrom].seq[start-1:end]) for start, end in coords)
    
    # Reverse complement if on negative strand
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    
    return seq

gencode_cds_coordinates['cds_sequence'] = gencode_cds_coordinates.apply(get_cds_sequence, axis=1)

```
Protein sequences were translated from the coding sequences of each transcript as:

```python
def cds_to_protein(cds_seq):
    if cds_seq is None or len(cds_seq) == 0:
        return None
    seq_obj = Seq(cds_seq)
    try:
        protein = str(seq_obj.translate(table=1, cds=True))
    except Exception as e:
        protein = str(seq_obj.translate(to_stop=True))
    return protein

gencode_cds_coordinates['protein_sequence'] = gencode_cds_coordinates['cds_sequence'].apply(cds_to_protein)
```

This generates the final **gencode_cds_coordinates** dataframe containing the protein sequences of all human coding sequences mapped to genomic coordinates. 

# Mapping protein sequences to genomic CDS coordinates

After generating the **gencode_cds_coordinates** dataframe, the **pro2gene** function is defined as: 

```python

def pro2gene(row, query_protein):
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

This function facilitates the use of a protein sequence query to identify its corresponding genomic coordinates.

An example is shown below:

```python

query_protein = "DLVILLYETALLSSGFSLEDPQTHANRIYRMIKLGLGIDEDDPTADDTSAAVTEEMPPLE"

gencode_cds_coordinates["matched_cds_coords"] = gencode_cds_coordinates.apply(
    lambda row: pro2gene(row, query_protein),
    axis=1
)

matches = gencode_cds_coordinates[gencode_cds_coordinates["matched_cds_coords"].notnull()]

matches[["gene_name", "seqname", "strand", "gene_id", "transcript_id", "protein_id", "matched_cds_coords"]]

```

Output:

```
	gene_name	seqname	strand	gene_id	transcript_id	protein_id	matched_cds_coords
1387198	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000216281.13	ENSP00000216281.8	[[[102081751, 102081821], [102082111, 10208221...
1387224	HSP90AA1	chr14	-	ENSG00000080824.19	ENST00000334701.11	ENSP00000335153.7	[[[102081751, 102081821], [102082111, 10208221..
```

We can check if the resulting coding sequences code for a given query peptide using the following function:

```python
def reconstruct_cds_from_coords(genome_seq_dict, chrom, cds_coords, strand):
    """
    genome_seq_dict: dict of chromosome sequences (SeqRecords)
    chrom: chromosome name
    cds_coords: list of [start, end] genomic coordinates
    strand: '+' or '-'
    """
    cds_seq = ""
    for start, end in cds_coords:
        # Access the .seq of SeqRecord and slice as a string
        cds_seq += str(genome_seq_dict[chrom].seq[start-1:end])
    if strand == '-':
        cds_seq = str(Seq(cds_seq).reverse_complement())
    return cds_seq


# Function to check if peptide exists in any translation frame
def peptide_in_cds_any_frame(cds_seq, query_protein):
    seq = Seq(cds_seq)
    frames = []

    # Forward frames
    for i in range(3):
        frames.append(str(seq[i:].translate(to_stop=False)))

    # Reverse complement frames
    rev_seq = seq.reverse_complement()
    for i in range(3):
        frames.append(str(rev_seq[i:].translate(to_stop=False)))

    # Check for peptide
    for f in frames:
        if query_protein in f:
            return True
    return False

# Example for the first row 
row = matches.iloc[0]
chrom = row['refseq_id']
strand = row['strand']
matched_coords_list = row['matched_cds_coords']  

# Loop over all matched segments 
found = False
for matched_coords in matched_coords_list:
    cds_seq = reconstruct_cds_from_coords(genome, chrom, matched_coords, strand)
    if peptide_in_cds_any_frame(cds_seq, query_protein):
        found = True
        break

print("Query protein found in matched CDS coordinates:", found)
```

Output:

```
Query protein found in matched CDS coordinates: True
```
