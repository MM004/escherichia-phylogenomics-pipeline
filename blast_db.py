import os

def create_blast_dbs(paths:dict):
    """ 2. Make blast databases from the sequences using makeblastdb. """

    orf_files_list = os.listdir(paths['orf_output_dir'])

    for orf_file in orf_files_list:
        if orf_file.endswith('_genes.fasta'):
            orf_path = os.path.join(paths['orf_output_dir'], orf_file)
            orf_prefix = orf_file.split('.fasta')[0]
            db_n_output_path = os.path.join(paths['blastn_db_dir'], orf_prefix)

            command2 = f'makeblastdb -in {orf_path} -dbtype nucl -out {db_n_output_path}'

            if not os.path.exists(db_n_output_path + '.nhr'):
                print('Running command:')
                print(command2)
                os.system(command2)

    ## proteins
    amino_acid_files_list = os.listdir(paths['amino_acid_output_dir'])

    for amino_acid_file in amino_acid_files_list:
        if amino_acid_file.endswith('_proteins.fasta'):
            amino_acid_path = os.path.join(paths['amino_acid_output_dir'], amino_acid_file)
            amino_acid_prefix = amino_acid_file.split('.fasta')[0]
            db_p_output_path = os.path.join(paths['blastp_db_dir'], amino_acid_prefix)

            command3 = f'makeblastdb -in {amino_acid_path} -dbtype prot -out {db_p_output_path}'

            if not os.path.exists(db_p_output_path + '.phr'):
                print('Running command:')
                print(command3)
                os.system(command3)

if __name__ == '__main__':
    print('Not intended to be run as a script! Import this module in another script!')