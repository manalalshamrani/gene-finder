from Bio import SeqIO
from Bio.Seq import Seq

def find_genes_in_fasta(file_path):
    # Parse the FASTA file
    for record in SeqIO.parse("GCA_000005845.2_ASM584v2_genomic.fna", "fasta"):
        sequence = str(record.seq)
        
        # Store found genes
        genes = []
        # Step 1: generate genes
        # Check each of the three reading frames
        for frame in range(3):
            # Read the sequence in the specified frame
            frame_sequence = sequence[frame:]
            
            # Find start codons
            start_index = frame_sequence.find("ATG")
            while start_index != -1:
                # Find stop codons after the start codon
                stop_index = -1
                for stop_codon in ["TAA", "TAG", "TGA"]:
                    temp_index = frame_sequence.find(stop_codon, start_index + 3)
                    if temp_index != -1:
                        if stop_index == -1 or temp_index < stop_index:
                            stop_index = temp_index
                
                # If a stop codon was found, extract the gene
                if stop_index != -1:
                    gene = frame_sequence[start_index:stop_index + 3]  # Include the stop codon
                    genes.append(gene)
                
                # Search for the next start codon
                start_index = frame_sequence.find("ATG", start_index + 1)
        
        # Print or return the genes found in the current record
        for gene in genes:
            print("Found genes in frame:",gene)

    # Step 2: Generate reverse complement of each gene found
	rev_comp_genes = [reverse_complement(gene) for gene in genes]
    
    # Step 3: Search for start and stop codons in the reverse complement
    	for rev_gene in rev_comp_genes:
		start_index = rev_gene.find("ATG")
		while start_index != -1:
			stop_index = -1
            # Find the closest stop codon
            for stop_codon in ["TAA", "TAG", "TGA"]:
                temp_index = rev_gene.find(stop_codon, start_index + 3)
                if temp_index != -1 and (stop_index == -1 or temp_index < stop_index):
                    stop_index = temp_index
            
            if stop_index != -1:
                gene = rev_gene[start_index:stop_index + 3]
                # Store or print the found gene from reverse complement
                print("Found gene in reverse complement:", gene)
            
            # Look for the next start codon
            start_index = rev_gene.find("ATG", start_index + 1)

    return genes, rev_comp_genes

# Function to generate reverse complement
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))


# Example usage
genes, rev_comp = find_genes_in_fasta("GCA_000005845.2_ASM584v2_genomic.fna")
