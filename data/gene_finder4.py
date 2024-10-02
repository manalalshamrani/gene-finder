from argparse import ArgumentParser
from Bio import SeqIO, Seq


def function_1(file_path: str):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences


def function_2_v1(sequence: str) \
        -> list:
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
    forward_orfs = function_2_v1(sequence=sequence)
    reverse_orfs = function_2_v1(sequence=sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])

    return forward_orfs + reverse_orfs


def function_3(orfs: list) \
        -> list:
    return list(set([Seq.Seq(orf[2]).translate(to_stop=True) for orf in orfs]))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    args = parser.parse_args()

    # Process the input file using function_1
    seqs = function_1(file_path=args.file)
    
    # Define the path to your output file
    output_file_path = "/home/alshammm/output.q4.txt"
    
    for key, value in seqs.items():
        # Find ORFs using function_2_v2 and filter by minimum length with function_3
        results = function_3(function_2_v2(sequence=value))
        
        # Open the output file once to write the results of this sequence
        with open(output_file_path, "a+") as file:
            # Write the number of ORFs found
            file.write("find " + str(len(results)) + " ORFs in sequence " + key + ":\n")
            
            # Write each result (ORF sequence) to the output file
            for result in results:
                file.write(str(result) + "\n")
