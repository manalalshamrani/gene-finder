import argparse
from Bio.Seq import Seq
import os

def read_sequence_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    sequence_data = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence_data

def find_open_reading_frames(sequence, ribosome_binding_site, scan_distance):
    stop_codons = ['TAA', 'TAG', 'TGA']
    reading_frames = []
    
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == 'ATG': 
                start = i
                upstream_region_start = max(0, start - scan_distance)
                upstream_region = sequence[upstream_region_start:start]
                if ribosome_binding_site in upstream_region:
                    for j in range(i + 3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:
                            reading_frames.append(sequence[start:j+3])
                            break
            i += 3
    return reading_frames

def find_all_reading_frames(sequence, ribosome_binding_site, scan_distance):
    reading_frames = []
    reading_frames.extend(find_open_reading_frames(sequence, ribosome_binding_site, scan_distance)) 
    
    reverse_complement_seq = str(Seq(sequence).reverse_complement())
    reading_frames.extend(find_open_reading_frames(reverse_complement_seq, ribosome_binding_site, scan_distance))  
    
    return reading_frames

def translate_reading_frames_to_proteins(reading_frames, min_protein_length):
    proteins = set()
    for reading_frame in reading_frames:
        protein = str(Seq(reading_frame).translate(to_stop=True))
        if len(protein) >= min_protein_length:
            proteins.add((protein, len(protein)))
    return proteins

def process_sequence_file(file_path, min_protein_length, ribosome_binding_site, scan_distance):
    sequence_data = read_sequence_file(file_path)
    reading_frames = find_all_reading_frames(sequence_data, ribosome_binding_site, scan_distance)
    proteins = translate_reading_frames_to_proteins(reading_frames, min_protein_length)
    
    file_dir = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    
    output_file = os.path.join(file_dir, f"{file_name}_output.txt")
    
    with open(output_file, 'w') as file:
        file.write(f"Proteins for file {file_path}:\n")
        for protein, length in proteins:
            file.write(f"{protein} {length}\n")
        file.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Extract ORFs, filter by length, RBS presence, and translate them into protein strings.")
    parser.add_argument('input_files', nargs='+', help="Paths to the input FASTA files containing the DNA sequences")
    parser.add_argument('--min_length', type=int, default=100, help="Minimum length of ORFs to be considered (in codons, default: 100)")
    parser.add_argument('--upstream_distance', type=int, default=20, help="Distance upstream of the start codon to search for RBS (default: 20)")
    parser.add_argument('--rbs_sequence', type=str, default='AGGAGG', help="Ribosome Binding Site (RBS) sequence to search for (default: AGGAGG)")
    
    args = parser.parse_args()
    
    for file_path in args.input_files:
        if "GCA" in file_path:
            process_sequence_file(file_path, args.min_length, args.rbs_sequence, args.upstream_distance)

if __name__ == '__main__':
    main()
