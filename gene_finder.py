from Bio.Seq import Seq
from Bio import SeqIO

def find_genes_in_fasta(file_path):
    # Parse the FASTA file
    for record in SeqIO.parse("GCA_000005845.2_ASM584v2_genomic.fna", "fasta"):
        sequence = str(record.seq)
        
        # Store found genes
        genes = []
        
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
        print(f"Found genes in frame {frame + 1}:")
        for gene in genes:
            print(gene)

# Example usage
find_genes_in_fasta("GCA_000005845.2_ASM584v2_genomic.fna")

