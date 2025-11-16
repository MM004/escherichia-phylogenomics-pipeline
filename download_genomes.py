from Bio import Entrez
import os

def download_g(paths:dir, genus = 'Escherichia', num_genomes = 10):
    """ 
    0. Download genomes from NCBI
    """
    # Set up your email for Entrez
    with open('email.txt', 'r') as e:
        Entrez.email = e.read()

    # Define the genus and the number of genomes to download
    genus = genus
    num_genomes = num_genomes  # Number of genomes to download

    # Search for bacterial genomes in NCBI
    search_query = f'{genus}[Organism] AND complete genome[Title]'

    with Entrez.esearch(db='nucleotide', term=search_query, retmax=num_genomes) as search_handle:
        search_results = Entrez.read(search_handle)
        genome_ids = search_results['IdList']

    if not genome_ids:
        print(f'No genomes found for genus {genus}.')
        exit()

    print(f'Found {len(genome_ids)} genomes.')

    # Download each genome in FASTA format
    download_dir = paths['input_genome_dir']
    for genome_id in genome_ids:
        with Entrez.efetch(db='nucleotide', id=genome_id, rettype='fasta', retmode='text') as fetch_handle:
            genome_data = fetch_handle.read()

        # Save the genome to a file
        genome_filename = os.path.join(download_dir, f'{genus}_{genome_id}.fasta')
        with open(genome_filename, 'w') as genome_file:
            genome_file.write(genome_data)

        print(f'Downloaded genome {genome_id} to {genome_filename}')

    print(f'All genomes downloaded and saved in {download_dir}.')

if __name__ == '__main__':
    print('Not intended to be run as a script! Import this module in another script!')