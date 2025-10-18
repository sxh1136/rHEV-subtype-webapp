import sys
import json
import os
import tempfile
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import subprocess
import numpy as np
import math


def calculate_k80_distance(fasta_file):
    # Read the aligned sequences
    alignment = AlignIO.read(fasta_file, "fasta")
    n = len(alignment)
    
    # Function to compute transitions and transversions between two sequences
    def transitions_transversions(seq1, seq2):
        transitions = 0
        transversions = 0
        transitions_pairs = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
        
        length = 0  # count valid comparisons (ignore gaps etc)
        for a, b in zip(seq1, seq2):
            if a in {'A','C','G','T'} and b in {'A','C','G','T'}:
                length += 1
                if a != b:
                    if (a,b) in transitions_pairs:
                        transitions += 1
                    else:
                        transversions += 1
        return transitions, transversions, length

    # Calculate distance matrix
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            ts, tv, length = transitions_transversions(alignment[i].seq, alignment[j].seq)
            if length == 0:
                dist = 0
            else:
                P = ts / length
                Q = tv / length
                # Kimura 2-parameter distance formula
                try:
                    dist = -0.5 * math.log(1 - 2*P - Q) - 0.25 * math.log(1 - 2*Q)
                except ValueError:
                    # if log argument <= 0 due to saturation, set distance to large number
                    dist = float('inf')
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    return dist_matrix

def main(input_fasta, existing_msa, output_dir):
    try:
        # Read input sequence
        input_seq = next(SeqIO.parse(input_fasta, "fasta"))

        # Generate a unique ID for the input sequence
        original_id = input_seq.id
        input_seq.id = f"QUERY_{original_id}"  # Add a prefix to ensure uniqueness
        input_seq.name = input_seq.id
        input_seq.description = input_seq.id
        new_input_id = input_seq.id

        # Save input sequence to a temporary file
        with tempfile.NamedTemporaryFile(mode='w+t', suffix=".fasta", delete=False) as temp_file:
            SeqIO.write([input_seq], temp_file.name, "fasta")
            new_seq_file = temp_file.name
        
        # Perform profile alignment using MUSCLE
        with tempfile.TemporaryDirectory() as tmp_dir:
            combined_alignment_file = os.path.join(tmp_dir, "combined.afa")
            muscle_cline = ["./muscle", "-profile", "-in1", existing_msa, "-in2", new_seq_file, "-out", combined_alignment_file]
            subprocess.run(muscle_cline, check=True, capture_output=True, text=True)
            
            # Read the combined alignment
            aligned_sequences = list(SeqIO.parse(combined_alignment_file, "fasta"))
        
            # Create a temporary file for the combined alignment to pass to R
            temp_alignment_file = os.path.join(tmp_dir, "alignment.fasta")
            SeqIO.write(aligned_sequences, temp_alignment_file, "fasta")

            # Calculate distance matrix using R function
            distance_matrix = calculate_k80_distance(temp_alignment_file)

        # Create a mapping of sequence IDs to their indices
        seq_ids = [seq.id for seq in aligned_sequences]
        id_to_index = {seq_id: idx for idx, seq_id in enumerate(seq_ids)}

        # Create a dictionary to store distances
        p_distances = {}
        for ref_seq in aligned_sequences:
            if ref_seq.id != new_input_id:
                # Get indices for the new input ID and reference ID
                new_index = id_to_index[new_input_id]
                ref_index = id_to_index[ref_seq.id]
                distance = distance_matrix[new_index, ref_index]
                p_distances[ref_seq.id] = distance

        # Find closest reference
        min_id, min_distance = min(p_distances.items(), key=lambda x: x[1])
        
        # Prepare output
        output = {
            "closest_reference": min_id,
            "p_distance": float(min_distance),  # Convert to native float
            "below_cutoff": bool(min_distance <= 0.1769),  # Convert to native bool
            "original_input_id": original_id  # Include the original ID in the output
        }

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Write JSON output
        with open(os.path.join(output_dir, "p_distance_output.json"), "w") as f:
            json.dump(output, f, indent=4)  # added indent for readability
            
        # Write distances file
        with open(os.path.join(output_dir, "p_distances.txt"), "w") as f:
            for ref_id, distance in p_distances.items():
                f.write(f"{ref_id}: {distance}\n")
                
        # Save the combined alignment to a file
        output_alignment_path = os.path.join(output_dir, "updated_alignment.fasta")
        with open(output_alignment_path, "w") as f:
            SeqIO.write(aligned_sequences, output_alignment_path, "fasta")
        
        return output_alignment_path  # returning path
    
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running MUSCLE: {e.stderr}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python msa_pdistance.py <input_fasta> <existing_msa> <output_dir>")
        sys.exit(1)
        
    alignment_path = main(sys.argv[1], sys.argv[2], sys.argv[3])
