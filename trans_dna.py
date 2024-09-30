from Bio import SeqIO, Seq
import sys

def read_fasta_file(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def find_orfs(sequence):
    start_codon, stop_codons = "ATG", ["TAA", "TAG", "TGA"]
    orfs = []
    for frame in range(3):
        start_positions = []
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                start_positions.append(i)
            elif codon in stop_codons:
                while start_positions:
                    start_pos = start_positions.pop(0)
                    orf = sequence[start_pos:i+3]
                    orfs.append(orf)
    return orfs

def find_reverse_orfs(sequence):
    reverse_sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return find_orfs(reverse_sequence)

def find_all_orfs(sequence):
    forward_orfs = find_orfs(sequence)
    reverse_orfs = find_reverse_orfs(sequence)
    return forward_orfs + reverse_orfs

def translate_orfs(orfs):
    translated_orfs = set()
    for orf in orfs:
        translated_orf = str(Seq.Seq(orf).translate(to_stop=True))
        translated_orfs.add(translated_orf)
    return list(translated_orfs)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python orf_finder.py <input_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    sequences = read_fasta_file(file_path)
    for key, value in sequences.items():
        orfs = find_all_orfs(value)
        translated_orfs = translate_orfs(orfs)
        for orf in translated_orfs:
            print(orf)
