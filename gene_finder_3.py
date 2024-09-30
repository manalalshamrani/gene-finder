from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def find_genes(sequence):
    genes = []
    # Check each of the three reading frames
    for frame in range(3):
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

    return genes

def find_genes_in_fasta(file_path):
    # Store all found genes
    all_genes = []

    # Parse the FASTA file
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)

        # Step 1: Find genes in the original sequence
        genes = find_genes(sequence)
        all_genes.extend(genes)

        # Step 2: Generate reverse complement of the sequence
        rev_comp_sequence = reverse_complement(sequence)

        # Step 3: Find genes in the reverse complement
        rev_genes = find_genes(rev_comp_sequence)
        all_genes.extend(rev_genes)

    return all_genes

# Example usage
genes = find_genes_in_fasta("GCA_000005845.2_ASM584v2_genomic.fna")
for gene in genes:
    print("Found gene:", gene)

