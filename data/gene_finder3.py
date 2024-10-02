from argparse import ArgumentParser
from Bio import SeqIO, Seq


def function_1():
    # question 1.
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    args = parser.parse_args()

    sequences = {}
    for record in SeqIO.parse(args.file, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences


def function_2_v1(sequence: str) \
        -> list:
    # question 1.
    start_codon, stop_codons, orfs = "ATG", ["TAA", "TAG", "TGA"], []

    # check three different ORFs
    for frame in range(3):
        # check the current ORF.
        start_positions = []
        for location in range(frame, len(sequence), 3):
            codon = str(sequence[location: location + 3])
            if codon == start_codon:
                start_positions.append(location)
            elif codon in stop_codons:
                while start_positions:
                    start_pos = start_positions.pop(0)
                    orf = sequence[start_pos:location + 3]
                    orfs.append((start_pos, location + 3, orf))

    return orfs


def function_2_v2(sequence: str) \
        -> list:
    # question 2.
    forward_orfs = function_2_v1(sequence=sequence)
    reverse_orfs = function_2_v1(sequence=sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])

    return forward_orfs + reverse_orfs


def function_3(orfs: list) \
        -> list:
    # question 3.
    return list(set([Seq.Seq(orf[2]).translate(to_stop=True) for orf in orfs]))


if __name__ == "__main__":
    seqs = function_1()
    for key, value in seqs.items():
        for result in function_3(function_2_v2(sequence=value)):
            print(result)
