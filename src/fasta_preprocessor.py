import os

# data cleaning, transformation, splitting, etc.
class FastaPreprocessor:
    def __init__(self, gene_name, star_allele_fasta_file_path=None, chromosome_fasta_file_path=None):
        self.gene_name = gene_name
        self.fasta_file_path = fasta_file_path
        self.chromosome_fasta_file_path = chromosome_fasta_file_path
        self.dict_haplo = {}


    # returns a dictionary of star alleles and the gene body
    def create_fasta_dict_from_star_alleles(self):
        current_key = None
        current_sequence = ""
        
        with open(self.fasta_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'PharmVar' in line:
                    print(line)
                    current_key = line.split()[0][1:]
                    current_sequence = ""
                elif current_key:
                    current_sequence += line
                    self.dict_haplo[current_key] = current_sequence
        
        return dict_haplo


    def create_context_seq(self, end_pos, seq_len, end_padding):
        real_end_pos = end_pos + end_padding
        #start_pos = end_pos - seq_len
        start_pos = end_padding - seq_len
        # Initialize variables
        current_pos = 1  # 1-based position
        extracted_sequence = ''
        
        # Open the gzipped file and read the sequence
        with gzip.open(self.chromosome_fasta_file_path, 'rt') as f:
            # Skip the header line
            f.readline()
            
            for line in f:
                line = line.strip()
                line_length = len(line)
                
                # Check if the current line contains the start position
                if current_pos + line_length >= start_pos:
                    start_index = start_pos - current_pos
                    end_index = min(real_end_pos - current_pos, line_length)
                    extracted_sequence += line[start_index:end_index]
                    
                    # Update the start position for the next iteration
                    start_pos = current_pos + line_length + 1
                
                # Update the current position
                current_pos += line_length
                
                # Break if the end position is reached
                if current_pos >= real_end_pos:
                    break
        
        # Save the extracted sequence to a file
        data_dir_name = f'/oak/stanford/groups/rbaltman/alptartici/hyenaDNA/hyena-dna/data/raw/{self.gene_name}'

        if not os.path.exists(data_dir_name):
            os.makedirs(data_dir_name)
        
        with open(f'{data_dir_name}/ref_context_{seq_len//1000}kb_{self.gene_name}.txt', 'w') as f:
            f.write(extracted_sequence)


    # converts variant description string to a dictionary
    def _generate_variants_dict_for_row(self, variation):
        """
        Generate a dictionary of variants based on a single 'variation' value.
        
        Parameters:
        - variation (str): The 'variation' value from a single row
        
        Returns:
        - dict: A dictionary containing all possible combinations of references and variants
        """
        variants_dict = {}
        
        # Step 1: Split the string by the '>' sign
        ref_var = variation.split('>')
        
        # Step 2: Mark the string to the left of '>' as the reference and split by ','
        references = ref_var[0].replace('-', '').split(',')
        
        # Step 3: Mark the string to the right of '>' as the variant and split by ','
        variants = ref_var[1].replace('-', '').split(',')
        
        # Step 4: Generate a dictionary with all possible combinations
        for ref in references:
            if ref not in variants_dict:
                variants_dict[ref] = []
            for var in variants:
                if var not in variants_dict[ref]:
                    variants_dict[ref].append(var)
                    
        # Step 5: Return the dictionary
        return variants_dict

    


    
        
        


    


        