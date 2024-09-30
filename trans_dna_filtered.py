from Bio import SeqIO, Seq
import argparse

def load_fasta_file(file_path):
    sequence_data = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_data[record.id] = str(record.seq)
    return sequence_data

def identify_open_reading_frames(sequence):
    start_codon, stop_codons = "ATG", ["TAA", "TAG", "TGA"]
    open_reading_frames = []
    for frame in range(3):
        start_positions = []
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                start_positions.append(i)
            elif codon in stop_codons:
                while start_positions:
                    start_pos = start_positions.pop(0)
                    open_reading_frame = sequence[start_pos:i+3]
                    open_reading_frames.append(open_reading_frame)
    return open_reading_frames

def identify_reverse_open_reading_frames(sequence):
    reverse_sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return identify_open_reading_frames(reverse_sequence)

def identify_all_open_reading_frames(sequence):
    forward_open_reading_frames = identify_open_reading_frames(sequence)
    reverse_open_reading_frames = identify_reverse_open_reading_frames(sequence)
    return forward_open_reading_frames + reverse_open_reading_frames

def filter_open_reading_frames_by_length(open_reading_frames, min_length):
    filtered_open_reading_frames = []
    for open_reading_frame in open_reading_frames:
        if len(open_reading_frame) >= min_length * 3:
            filtered_open_reading_frames.append(open_reading_frame)
    return filtered_open_reading_frames

def translate_open_reading_frames(open_reading_frames):
    translated_open_reading_frames = set()
    for open_reading_frame in open_reading_frames:
        translated_open_reading_frame = str(Seq.Seq(open_reading_frame).translate(to_stop=True))
        translated_open_reading_frames.add(translated_open_reading_frame)
    return list(translated_open_reading_frames)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find Open Reading Frames in a DNA sequence")
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    parser.add_argument("-l", "--length", help="minimum length of Open Reading Frame", required=True)
    args = parser.parse_args()

    sequence_data = load_fasta_file(args.file)
    for key, value in sequence_data.items():
        open_reading_frames = identify_all_open_reading_frames(value)
        filtered_open_reading_frames = filter_open_reading_frames_by_length(open_reading_frames, int(args.length))
        translated_open_reading_frames = translate_open_reading_frames(filtered_open_reading_frames)
        for translated_open_reading_frame in translated_open_reading_frames:
            print(translated_open_reading_frame)
