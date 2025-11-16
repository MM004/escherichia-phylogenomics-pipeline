import os
from Bio import SeqIO

def muscle_alignment(paths:dir):
    """ 5. For each of the fasta files of homologous sequences, run MUSCLE to align them. """
    
    # nucleotides
    ## contactenate all the homologous genes
    with open(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta'), 'w') as concatenated_handle:
        for fasta_file in os.listdir(paths['homologous_genes_dir']):
            if fasta_file.endswith('.fasta'):
                fasta_path = os.path.join(paths['homologous_genes_dir'], fasta_file)
                with open(fasta_path, 'r') as fasta_handle:
                    concatenated_handle.write(fasta_handle.read())

    ## check for duplicated sequences in concatenated_genes.fasta and keep only one
    with open(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta'), 'r') as concatenated_handle:
        records = list(SeqIO.parse(concatenated_handle, 'fasta'))
        # Create a dictionary to retain only the first occurrence of each record ID
        unique_records = {}
        for record in records:
            if record.id not in unique_records:
                unique_records[record.id] = record  # Retain the first occurrence
        # Write the deduplicated records to a new file
        with open(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes_new.fasta'), 'w') as output_handle:
            SeqIO.write(unique_records.values(), output_handle, 'fasta')

    # Delete the old concatenated file
    if os.path.exists(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta')):
        os.remove(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta'))
        print(f"Old concatenated file '{os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta')}' has been deleted.")
    else:
        print(f"File '{os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta')}' not found. Skipping deletion.")

    # Rename the new concatenated file
    os.rename(os.path.join(paths['homologous_genes_dir'], 'concatenated_genes_new.fasta'), os.path.join(paths['homologous_genes_dir'], 'concatenated_genes.fasta'))

    ## allign all the homologous genes and the concatenated genes
    for fasta_file in os.listdir(paths['homologous_genes_dir']):
        if fasta_file.endswith('.fasta'):
            fasta_path = os.path.join(paths['homologous_genes_dir'], fasta_file)
            muscle_output_path = os.path.join(paths['aligned_genes_dir'], fasta_file).replace(".fasta", ".afa")
            command6 = f'muscle -align {fasta_path} -output {muscle_output_path}'

            if not os.path.exists(muscle_output_path):
                print('Running command:')
                print(command6)
                os.system(command6)
    ###############################################################################################################
    #  
    # proteins
    ## contactenate all the homologous proteins
    with open(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta'), 'w') as concatenated_handle:
        for fasta_file in os.listdir(paths['homologous_proteins_dir']):
            if fasta_file.endswith('.fasta'):
                fasta_path = os.path.join(paths['homologous_proteins_dir'], fasta_file)
                with open(fasta_path, 'r') as fasta_handle:
                    concatenated_handle.write(fasta_handle.read())

    ## check for duplicated sequences in concatenated_proteins.fasta and keep only one
    with open(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta'), 'r') as concatenated_handle:
        records = list(SeqIO.parse(concatenated_handle, 'fasta'))
        # Create a dictionary to retain only the first occurrence of each record ID
        unique_records = {}
        for record in records:
            if record.id not in unique_records:
                unique_records[record.id] = record  # Retain the first occurrence
        # Write the deduplicated records to a new file
        with open(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins_new.fasta'), 'w') as output_handle:
            SeqIO.write(unique_records.values(), output_handle, 'fasta')

    # Delete the old concatenated file
    if os.path.exists(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta')):
        os.remove(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta'))
        print(f"Old concatenated file '{os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta')}' has been deleted.")
    else:
        print(f"File '{os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta')}' not found. Skipping deletion.")
    
    # Rename the new concatenated file
    os.rename(os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins_new.fasta'), os.path.join(paths['homologous_proteins_dir'], 'concatenated_proteins.fasta'))
    
    ## allign all the homologous proteins and the concatenated proteins
    for fasta_file in os.listdir(paths['homologous_proteins_dir']):
        if fasta_file.endswith('.fasta'):
            fasta_path = os.path.join(paths['homologous_proteins_dir'], fasta_file)
            muscle_output_path = os.path.join(paths['aligned_proteins_dir'], fasta_file).replace(".fasta", ".afa")
            command7 = f'muscle -align {fasta_path} -output {muscle_output_path}'

            if not os.path.exists(muscle_output_path):
                print('Running command:')
                print(command7)
                os.system(command7)