import os

def prodigal_prediction(paths:dict):
    """ 1. Predict genes using prodigal and export the predicted amino acid and gene 
        (open reading frame at the nucleotide level) sequences using prodigal.  
    """
    ## extension of genome files
    genome_ext = '.fasta'
    genomes_list = os.listdir(paths['input_genome_dir'])

    for genome in genomes_list:
        if genome.endswith(genome_ext):
            genome_path = os.path.join(paths['input_genome_dir'], genome)
            genome_prefix = genome.split(genome_ext)[0]
            orf_output_path = os.path.join(paths['orf_output_dir'], genome_prefix + '_genes.fasta')
            amino_acid_output_path = os.path.join(paths['amino_acid_output_dir'], genome_prefix + '_proteins.fasta')

            command1 = f'prodigal -i {genome_path} -a {amino_acid_output_path} -d {orf_output_path}'

            if not os.path.exists(amino_acid_output_path) or not os.path.exists(orf_output_path):
                print('Running command:')
                print(command1)
                os.system(command1)
                
if __name__ == '__main__':
    print('Not intended to be run as a script! Import this module in another script!')