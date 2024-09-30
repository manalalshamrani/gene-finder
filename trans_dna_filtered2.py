import sys

# Default RBS sequences and search parameters
RBS_SEQUENCES = ["AGGAGG"]  # Add more RBS sequences as needed
MIN_CODON_LENGTH = 100
SEARCH_WINDOW = 20

def find_codons(genome_sequence):
    codons = []
    for i in range(len(genome_sequence) - 2):
        if genome_sequence[i:i+3] in ["ATG", "GTG", "TTG"]:  # Start codons
            for j in range(i+3, len(genome_sequence), 3):
                if genome_sequence[j:j+3] in ["TAA", "TAG", "TGA"]:  # Stop codons
                    codons.append((i, j+3))
                    break
    return codons

def filter_codons_by_binding(codons, genome_sequence):
    valid_codons = []
    for start, end in codons:
        orf = genome_sequence[start:end]
        if len(orf) % 3 == 0 and len(orf) >= 300:
            valid_codons.append((start, end))
    return valid_codons

def run():
    genome_sequence = ""
    with open(sys.argv[1], "r") as file:
        for line in file:
            if not line.startswith(">"):
                genome_sequence += line.strip()

    codons = find_codons(genome_sequence)
    valid_codons = filter_codons_by_binding(codons, genome_sequence)

    if valid_codons:
        for start, end in valid_codons:
            orf = genome_sequence[start:end]
            print(f"Valid ORF: Start: {start}, End: {end}, Sequence: {orf}")
    else:
        print("No valid ORFs found.")

if __name__ == "__main__":
    import sys
    run()
