# Escherichia Phylogenomics Pipeline

**Author:** Mathias Mayrgündter  

A Python-based workflow that downloads bacterial genomes, predicts genes, identifies homologs, aligns them, and constructs phylogenetic trees at both the protein and nucleotide levels.

---

## Pipeline Summary

This script performs all steps required for the phylogenomic assignment:

- Download genomes from NCBI for a chosen genus  
- Predict genes using Prodigal (protein sequences + nucleotide ORFs)  
- Create BLAST databases using `makeblastdb`  
- Run all-vs-all BLAST (`blastp` and `blastn`)  
- Parse BLAST tables and automatically select homologous protein and gene sets  
- Align sequences with MUSCLE  
- Build phylogenetic trees with IQ-TREE (individual and concatenated)  
- Output directories are automatically created under `dataflow/`

---

## Dependencies

Install via conda (recommended):

```bash
conda install -c bioconda prodigal blast muscle iqtree
pip install biopython pandas tqdm
```

---

## Running the Pipeline

```bash
python run_pipeline.py
```

**Default settings:**
- Genus: *Escherichia*  
- 10 genomes  
- Output folder: `dataflow/`  

Edit variables at the top of the script if you want to change genus or genome count.

---

## Output Structure

```
dataflow/
    01-genomes/
    nucleotides/
        02-orfs/
        03-blastn-db/
        04-blastn-out/
        05-homologous_genes/
        06-aligned_genes/
        07-phylogenetic_trees/
    proteins/
        02-amino_acids/
        03-blastp-db/
        04-blastp-out/
        05-homologous_proteins/
        06-aligned_proteins/
        07-phylogenetic_trees/
```

---

## Notes

- You must visualize the trees and run BLAST (online) to annotate representative sequences  
- Example tree images + short interpretation see pdf Documents  
