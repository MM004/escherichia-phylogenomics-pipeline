import os
import subprocess

def blast_all_vs_all(paths:dict):
    """ 3. For the following step, you will compare all the genomes against each other (‘all vs all’) 
    using blastp and blastn to identify sets of homologous proteins and genes (open reading frame at the nucleotide level) """

    orf_files_list = os.listdir(paths['orf_output_dir'])
    amino_acid_files_list = os.listdir(paths['amino_acid_output_dir'])

    # nucleotides
    ## create log file for blastn
    with open(os.path.join(paths['blastn_output_dir'], 'blastn.log'), 'w') as log_file:
        log_file.write('Blastn log file\n')
    # iterate over all orf files
    for orf_file in orf_files_list:
        orf_prefix = orf_file.split('.fasta')[0]
        orf_file_path = os.path.join(paths['orf_output_dir'], orf_file)
        print(f'Running blastn for {orf_prefix}')
        print('='*50)

        for blastn_db_file in os.listdir(paths['blastn_db_dir']):
            if blastn_db_file.endswith('.nhr'):
                blastn_db_prefix = blastn_db_file.split('.nhr')[0]
                blastn_db_path = os.path.join(paths['blastn_db_dir'], blastn_db_prefix)
                blastn_output_path = os.path.join(paths['blastn_output_dir'], orf_prefix.split('_')[1] + '_VS_' +  blastn_db_prefix.split('_')[1] + '.txt')

                command4 = f'blastn -query {orf_file_path} -db {blastn_db_path} -max_target_seqs 1 -evalue 1E-5 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out {blastn_output_path}'

                if not os.path.exists(blastn_output_path):
                    with open(os.devnull, 'w') as devnull:  # suppress output
                        print('Comparing', orf_prefix, 'against', blastn_db_prefix)
                        subprocess.run(command4, shell=True, stdout=devnull, stderr=subprocess.STDOUT)

        ## append to log file that blastn has finished running for this orf file
        with open(os.path.join(paths['blastn_output_dir'], 'blastn.log'), 'a') as log_file:
            log_file.write(f"Finished running blastn for {orf_prefix}\n")
        ## print that blastn has finished running for this orf file
        print('-'*50)
        print(f'Finished running blastn for {orf_prefix}')

    ## finish blastn
    with open(os.path.join(paths['blastn_output_dir'], 'blastn.log'), 'a') as log_file:
            log_file.write(f"Finished running blastn for all orf files\n")

    # proteins
    ## create log file for blastp
    with open(os.path.join(paths['blastp_output_dir'], 'blastp.log'), 'w') as log_file:
        log_file.write('Blastp log file\n')
    ## iterate over all amino acid files
    for amino_acid_file in amino_acid_files_list:
        amino_acid_file_path = os.path.join(paths['amino_acid_output_dir'], amino_acid_file)
        amino_acid_prefix = amino_acid_file.split('.fasta')[0]
        print(f'Running blastp for {amino_acid_prefix}')
        print('='*50)

        for blastp_db_file in os.listdir(paths['blastp_db_dir']):
            if blastp_db_file.endswith('.phr'):
                blastp_db_prefix = blastp_db_file.split('.phr')[0]
                blastp_db_path = os.path.join(paths['blastp_db_dir'], blastp_db_prefix)
                blastp_output_path = os.path.join(paths['blastp_output_dir'], amino_acid_prefix.split('_')[1] + '_VS_' +  blastp_db_prefix.split('_')[1] + '.txt')

                command5 = f'blastp -query {amino_acid_file_path} -db {blastp_db_path} -max_target_seqs 1 -evalue 1E-5 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out {blastp_output_path}'

                if not os.path.exists(blastp_output_path):
                    with open(os.devnull, 'w') as devnull:  # suppress output
                        print('Comparing', amino_acid_prefix, 'against', blastp_db_prefix)
                        subprocess.run(command5, shell=True, stdout=devnull, stderr=subprocess.STDOUT)
                    os.system(command5)

        ## append to log file that blastp has finished running for this amino acid file
        with open(os.path.join(paths['blastp_output_dir'], 'blastp.log'), 'a') as log_file:
            log_file.write(f"Finished running blastp for {amino_acid_prefix}\n")
        ## print that blastp has finished running for this amino acid file
        print('-'*50)
        print(f'Finished running blastp for {amino_acid_prefix}')

    ## finish blastp
    with open(os.path.join(paths['blastp_output_dir'], 'blastp.log'), 'a') as log_file:
            log_file.write(f"Finished running blastp for all amino acid files\n")

if __name__ == '__main__':
    print('Not intended to be run as a script! Import this module in another script!')