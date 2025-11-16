from Bio import SeqIO
import pandas as pd
import os

def parse_blast_output(paths:dir):
    """  Parse the output tables to select homologous sequences based on cut-offs that you decide 
    are appropriate. State why you selected these cut-offs. 
    Select sets of sequences separately for homologs at the amino acid (5 sets) and nucleotide (5 sets) 
    level (this selection can be done randomly). Then, combine the sequences for each set into individual fasta files. 
    To compile the individual fasta files with homologous sequences, you may do this ‘by hand’ (i.e., copy and paste) and will receive 50% of the points for this step. 
    If you automate selecting and compiling the fasta files in the python script, you will receive full points.  """
    col_headers = ["qseqid", "sseqid", "pident","evalue", "sstart", "send", "length"]

    ## nucleotides
    output_parse_files_blastn = sorted([f for f in os.listdir(paths['blastn_output_dir']) if f.endswith('.txt')])
    list_df_tmp = []

    for output_file in output_parse_files_blastn:
        input_path = os.path.join(paths['blastn_output_dir'], output_file)
        
        # Read in the blast output like R: no quoting, strip whitespace from fields
        df_blast_tmp = pd.read_csv(
            input_path,
            sep="\t",
            header=None,
            names=col_headers,
            engine='python',  
            skipinitialspace=True  
        )
        # Filter the blast output to keep only the hits with pident >= 96 and pident < 100
        # because we want to select homologous genes and not identical genes
        df_blast_tmp = df_blast_tmp[df_blast_tmp['pident'] >= 96]
        df_blast_tmp = df_blast_tmp[df_blast_tmp['pident'] < 100]
        df_blast_tmp['file'] = output_file

        list_df_tmp.append(df_blast_tmp)

    df_blast_clean_compiled = pd.concat(list_df_tmp, ignore_index=True)

    # Separate file into query_genes_file and subject_genes_file
    df_blast_clean_metrics = df_blast_clean_compiled.copy()
    df_blast_clean_metrics[['query_genes_file', 'subject_genes_file']] = df_blast_clean_metrics['file'].str.split('_VS_', expand=True)
    df_blast_clean_metrics['subject_genes_file'] = df_blast_clean_metrics['subject_genes_file'].str.replace('.txt$', '', regex=True)

    # Calculate n_genomes_hit
    df_blast_clean_metrics['n_genomes_hit'] = df_blast_clean_metrics.groupby('qseqid')['subject_genes_file'].transform('nunique')

    # Calculate mean_pident
    df_blast_clean_metrics['mean_pident'] = df_blast_clean_metrics.groupby('qseqid')['pident'].transform('mean')

    # Select gene sets: n_genomes_hit == 9 and Bottom 5 by mean_pident
    df_blast_clean_select_gene_sets = (
        df_blast_clean_metrics[df_blast_clean_metrics['n_genomes_hit'] == 9]
        .loc[:, ['qseqid', 'mean_pident']]
        .drop_duplicates()
        .sort_values('mean_pident', ascending=True)
        .tail(5)
    )

    # Print out gene sets for each top gene and extract sequences for the homologous genes into fasta files
    for gene in df_blast_clean_select_gene_sets['qseqid'].unique():
        single_gene_set_blast_df = (
            df_blast_clean_metrics[df_blast_clean_metrics['qseqid'] == gene]
            .loc[:, ['qseqid', 'sseqid']]
            .drop_duplicates()
        )
        # You need to compile fasta files with these sets of homologous genes
        gene_set = single_gene_set_blast_df['sseqid'].unique().tolist()
        print(gene, gene_set)
        # Extract sequences for the homologous genes
        output_fasta_file = os.path.join(paths['homologous_genes_dir'], f'{gene}.fasta')
        if not os.path.exists(output_fasta_file):
            with open(output_fasta_file, 'w') as output_handle:
                for fasta_file in os.listdir(paths['orf_output_dir']):
                    if fasta_file.endswith(".fasta"):
                        fasta_path = os.path.join(paths['orf_output_dir'], fasta_file)
                        for record in SeqIO.parse(fasta_path, "fasta"):
                            if record.id in gene_set or record.id == gene:
                                SeqIO.write(record, output_handle, "fasta")
                                print(f"Extracted: {record.id}")
    
    ############################################################################################################
    output_parse_files_blastp = sorted([f for f in os.listdir(paths['blastp_output_dir']) if f.endswith('.txt')])
    list_df_tmp = []

    for output_file in output_parse_files_blastp:
        input_path = os.path.join(paths['blastp_output_dir'], output_file)
        
        # Read in the blast output like R: no quoting, strip whitespace from fields
        df_blast_tmp = pd.read_csv(
            input_path,
            sep="\t",
            header=None,
            names=col_headers,
            engine='python',  
            skipinitialspace=True  
        )
        # Filter the blast output to keep only the hits with pident >= 96 and pident < 100
        # because we want to select homologous proteins and not identical proteins
        df_blast_tmp = df_blast_tmp[df_blast_tmp['pident'] >= 96]
        df_blast_tmp = df_blast_tmp[df_blast_tmp['pident'] < 100]
        df_blast_tmp['file'] = output_file

        list_df_tmp.append(df_blast_tmp)

    df_blast_clean_compiled = pd.concat(list_df_tmp, ignore_index=True)

    # Separate file into query_proteins_file and subject_proteins_file
    df_blast_clean_metrics = df_blast_clean_compiled.copy()
    df_blast_clean_metrics[['query_proteins_file', 'subject_proteins_file']] = df_blast_clean_metrics['file'].str.split('_VS_', expand=True)
    df_blast_clean_metrics['subject_proteins_file'] = df_blast_clean_metrics['subject_proteins_file'].str.replace('.txt$', '', regex=True)

    # Calculate n_genomes_hit
    df_blast_clean_metrics['n_genomes_hit'] = df_blast_clean_metrics.groupby('qseqid')['subject_proteins_file'].transform('nunique')

    # Calculate mean_pident
    df_blast_clean_metrics['mean_pident'] = df_blast_clean_metrics.groupby('qseqid')['pident'].transform('mean')

    # Select protein sets: n_genomes_hit == 9 and Bottom 5 by mean_pident
    df_blast_clean_select_protein_sets = (
        df_blast_clean_metrics[df_blast_clean_metrics['n_genomes_hit'] == 9]
        .loc[:, ['qseqid', 'mean_pident']]
        .drop_duplicates()
        .sort_values('mean_pident', ascending=True)
        .tail(5)
    )

    # Print out protein sets for each top protein and extract sequences for the homologous proteins ino fasta files
    for protein in df_blast_clean_select_protein_sets['qseqid'].unique():
        single_protein_set_blast_df = (
            df_blast_clean_metrics[df_blast_clean_metrics['qseqid'] == protein]
            .loc[:, ['qseqid', 'sseqid']]
            .drop_duplicates()
        )
        # You need to compile fasta files with these sets of homologous proteins
        protein_set = single_protein_set_blast_df['sseqid'].unique().tolist()
        print(protein, protein_set)
        # Extract sequences for the homologous proteins
        output_fasta_file = os.path.join(paths['homologous_proteins_dir'], f'{protein}.fasta')
        if not os.path.exists(output_fasta_file):
            with open(output_fasta_file, 'w') as output_handle:
                for fasta_file in os.listdir(paths['amino_acid_output_dir']):
                    if fasta_file.endswith(".fasta"):
                        fasta_path = os.path.join(paths['amino_acid_output_dir'], fasta_file)
                        for record in SeqIO.parse(fasta_path, "fasta"):
                            if record.id in protein_set or record.id == protein:
                                SeqIO.write(record, output_handle, "fasta")
                                print(f"Extracted: {record.id}")