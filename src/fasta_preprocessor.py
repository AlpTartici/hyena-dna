import os



# data cleaning, transformation, splitting, etc
class FastaPreprocessor:
    def __init__(self, gene_name, fasta_file_path=None):
        self.gene_name = gene_name
        self.fasta_file_path = fasta_file_path


    # returns a dictionary of star alleles and the gene body
    def create_fasta_dict_from_star_alleles():
        dict_haplo = {}
        current_key = None
        current_sequence = ""
        
        file_path = '/oak/stanford/groups/rbaltman/alptartici/hyenaDNA/hyena-dna/alp_data/CYP2C9.haplotypes.fasta'
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'PharmVar' in line:
                    print(line)
                    current_key = line.split()[0][1:]
                    current_sequence = ""
                elif current_key:
                    current_sequence += line
                    dict_haplo[current_key] = current_sequence
        
        return dict_haplo


    


        