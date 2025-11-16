import os
from Bio import SeqIO

def iqtree_constructer(paths:dir):
    """ 6. Construct a phylogenetic tree using iqtree. This step can be conducted online to receive 
    50% of the points (http://www.iqtree.org/). If you automatically generate the trees with your 
    python script you will receive full points. You can find example commands for running 
    iqtree on their website. Carry out this step (6) with all individual alignments 
    and for a concatenated alignment of all protein and gene alignments.  """

    # nucleotides
    for aligned_file in os.listdir(paths['aligned_genes_dir']):
        if aligned_file.endswith('.afa'):
            aligned_path = os.path.join(paths['aligned_genes_dir'], aligned_file)
            iqtree_output_path = os.path.join(paths['iqtree_genes_dir'], aligned_file).replace(".afa", "")
            command8 = f'iqtree2 -s {aligned_path} -m MFP -bb 1000 -nt AUTO -st DNA -pre {iqtree_output_path}'

            if not os.path.exists(iqtree_output_path + ".treefile"):
                print('Running command:')
                print(command8)
                os.system(command8)

    ## check if the tree file exists otherwise run command without bootstrapping
    for aligned_file in os.listdir(paths['aligned_genes_dir']):
        aligned_path = os.path.join(paths['aligned_genes_dir'], aligned_file)
        iqtree_output_path = os.path.join(paths['iqtree_genes_dir'], aligned_file).replace(".afa", "")
        if aligned_file.endswith('.afa'):

            if not os.path.exists(iqtree_output_path + ".treefile"):

                command8a = f'iqtree2 -s {aligned_path} -m MFP -nt AUTO -st DNA -pre {iqtree_output_path}'

                print('Running command:')
                print(command8a)
                os.system(command8a)

    # proteins
    for aligned_file in os.listdir(paths['aligned_proteins_dir']):
        if aligned_file.endswith('.afa'):
            aligned_path = os.path.join(paths['aligned_proteins_dir'], aligned_file)
            iqtree_output_path = os.path.join(paths['iqtree_proteins_dir'], aligned_file).replace(".afa", "")
            command9 = f'iqtree2 -s {aligned_path} -m MFP -bb 1000 -nt AUTO -st AA -pre {iqtree_output_path}'

            if not os.path.exists(iqtree_output_path + ".treefile"):
                print('Running command:')
                print(command9)
                os.system(command9)


    ## check if the tree file exists otherwise run command without bootstrapping
    for aligned_file in os.listdir(paths['aligned_proteins_dir']):
        if aligned_file.endswith('.afa'):
            aligned_path = os.path.join(paths['aligned_proteins_dir'], aligned_file)
            iqtree_output_path = os.path.join(paths['iqtree_proteins_dir'], aligned_file).replace(".afa", "")

            if not os.path.exists(iqtree_output_path + ".treefile"):

                command9a = f'iqtree2 -s {aligned_path} -m MFP -nt AUTO -st AA -pre {iqtree_output_path}'

                print('Running command:')
                print(command9a)
                os.system(command9a)
                