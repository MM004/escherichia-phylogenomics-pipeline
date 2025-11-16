import os
from download_genomes import download_g
from prediction import prodigal_prediction
from blast_db import create_blast_dbs
from blast import blast_all_vs_all
from parse import parse_blast_output
from alignment import muscle_alignment
from phylogenetic_tree import iqtree_constructer

# =============================================================================
# =============================================================================
# Create a directory for dataflow
dataflow_dirs = ['dataflow', 'dataflow/nucleotides', 'dataflow/proteins',
                 'dataflow/01-genomes','dataflow/nucleotides/02-orfs','dataflow/proteins/02-amino_acids',
                 'dataflow/nucleotides/03-blastn-db', 'dataflow/proteins/03-blastp-db',
                 'dataflow/nucleotides/04-blastn-out', 'dataflow/proteins/04-blastp-out',
                 'dataflow/nucleotides/05-homologous_genes', 'dataflow/proteins/05-homologous_proteins',
                 'dataflow/nucleotides/06-aligned_genes', 'dataflow/proteins/06-aligned_proteins',
                 'dataflow/nucleotides/07-phylogenetic_trees', 'dataflow/proteins/07-phylogenetic_trees']

for dir in dataflow_dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)

# =============================================================================
# =============================================================================
# directories 
input_genome_dir = 'dataflow/01-genomes/'

## nucleotides
orf_output_dir = 'dataflow/nucleotides/02-orfs'
blastn_db_dir = 'dataflow/nucleotides/03-blastn-db'
blastn_output_dir = 'dataflow/nucleotides/04-blastn-out'
homologous_genes_dir = 'dataflow/nucleotides/05-homologous_genes'
aligned_genes_dir = 'dataflow/nucleotides/06-aligned_genes'
iqtree_genes_dir = 'dataflow/nucleotides/07-phylogenetic_trees'

## proteins
amino_acid_output_dir = 'dataflow/proteins/02-amino_acids'
blastp_db_dir = 'dataflow/proteins/03-blastp-db'
blastp_output_dir = 'dataflow/proteins/04-blastp-out'
homologous_proteins_dir = 'dataflow/proteins/05-homologous_proteins'
aligned_proteins_dir = 'dataflow/proteins/06-aligned_proteins'
iqtree_proteins_dir = 'dataflow/proteins/07-phylogenetic_trees'

## dictionary of paths
var_paths_dir = {'input_genome_dir': input_genome_dir, 
                    'orf_output_dir': orf_output_dir, 
                    'blastn_db_dir': blastn_db_dir,
                    'blastn_output_dir': blastn_output_dir,
                    'homologous_genes_dir': homologous_genes_dir,
                    'aligned_genes_dir': aligned_genes_dir,
                    'iqtree_genes_dir': iqtree_genes_dir,
                    'amino_acid_output_dir': amino_acid_output_dir,
                    'blastp_db_dir': blastp_db_dir,
                    'blastp_output_dir': blastp_output_dir,
                    'homologous_proteins_dir': homologous_proteins_dir,
                    'aligned_proteins_dir': aligned_proteins_dir,
                    'iqtree_proteins_dir': iqtree_proteins_dir}

# =============================================================================
# =============================================================================
# Modules
## 0. Download genomes from NCBI
download_g(var_paths_dir)

## 1. Predict genes using prodigal and export the predicted amino acid and gene 
## (open reading frame at the nucleotide level) sequences using prodigal. 
prodigal_prediction(var_paths_dir)

## 2. Make blast databases from the sequences using makeblastdb.
create_blast_dbs(var_paths_dir)

## 3. For the following step, you will compare all the genomes against each other (‘all vs all’) 
## using blastp and blastn to identify sets of homologous proteins and genes (open reading frame at the nucleotide level
blast_all_vs_all(var_paths_dir)

## 4. Parse the output tables to select homologous sequences based on cut-offs that you decide are appropriate. 
## State why you selected these cut-offs. Select sets of sequences separately for homologs at the amino acid (5 sets) and nucleotide (5 sets) level 
## (this selection can be done randomly). Then, combine the sequences for each set into individual fasta files.
##  To compile the individual fasta files with homologous sequences, you may do this ‘by hand’ (i.e., copy and paste) 
## and will receive 50% of the points for this step. If you automate selecting and compiling the fasta files in the python script, you will receive full points. 
parse_blast_output(var_paths_dir)

## 5. For each of the fasta files of homologous sequences, run MUSCLE to align them.
muscle_alignment(var_paths_dir)

## 6. Construct a phylogenetic tree using iqtree. This step can be conducted online to receive 
## 50% of the points (http://www.iqtree.org/). If you automatically generate the trees with your 
## python script you will receive full points. You can find example commands for running 
## iqtree on their website. Carry out this step (6) with all individual alignments 
## and for a concatenated alignment of all protein and gene alignments. 
iqtree_constructer(var_paths_dir)