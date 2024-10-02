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


def function_3_v1(orfs: list,
                  minimum_length: int) \
        -> list:
    sequences = set()
    for orf in orfs:
        sequence = str(Seq.Seq(orf[2]).translate(to_stop=True))
        if len(sequence) >= minimum_length:
            sequences.add(sequence)
    return list(sequences)


def function_3_v2(sequence: str,
                  rbs_list: list,
                  upstream: int,
                  minimum_length: int) \
        -> list:
    forward_sequence = sequence
    reverse_sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    sequences = set()
    for start, stop, frame in function_2_v1(sequence=forward_sequence):
        rbs_region, flag = forward_sequence[max(0, start - upstream): start], False
        for rbs in rbs_list:
            if rbs in rbs_region:
                flag = True
                break
        if flag:
            frame = str(Seq.Seq(frame).translate(to_stop=True))
            if len(frame) >= minimum_length:
                sequences.add(frame)

    for start, stop, frame in function_2_v1(sequence=reverse_sequence):
        rbs_region, flag = reverse_sequence[max(0, start - upstream): start], False
        for rbs in rbs_list:
            if rbs in rbs_region:
                flag = True
                break
        if flag:
            frame = str(Seq.Seq(frame).translate(to_stop=True))
            if len(frame) >= minimum_length:
                sequences.add(frame)

    return list(sequences)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    parser.add_argument("-l", "--length", help="minimum length of ORF", required=True)
    parser.add_argument("-u", "--upstream", help="the upstream distance to search for RBS", required=True)
    parser.add_argument("-r", "--rbs", help="RBS sequences", required=True)
    args = parser.parse_args()

    seqs = function_1(file_path=args.file)
    for key, value in seqs.items():
        results = function_3_v2(sequence=value,
                                rbs_list=args.rbs.split(","), upstream=int(args.upstream),
                                minimum_length=int(args.length))
        with open("/home/alshammm/output.q6" + args.file[48:61] + ".txt", "a+") as file:
            file.write("find " + str(len(results)) + " ORFs:" + "\n")
        for result in results:
            with open("/home/alshammm/output.q6" + args.file[48:61] + ".txt", "a+") as file:
                file.write(str(result) + "\n")
