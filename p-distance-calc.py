import sys
import json
import os
import tempfile
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
import subprocess

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
        
        # Create MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(aligned_sequences)
        
        # Calculate distance matrix 
        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(alignment)
        
        # Extract distances for input sequence
        p_distances = {}

        for seq in aligned_sequences:
            if seq.id != new_input_id:
                ref_id = seq.id
                distance = distance_matrix[new_input_id, ref_id]
                p_distances[ref_id] = distance
        
        # Find closest reference
        min_id, min_distance = min(p_distances.items(), key=lambda x: x[1])
        
        # Prepare output
        output = {
            "closest_reference": min_id,
            "p_distance": min_distance,
            "below_cutoff": min_distance <= 0.1833,
            "original_input_id": original_id  # Include the original ID in the output
        }

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Write JSON output
        with open(os.path.join(output_dir, "p_distance_output.json"), "w") as f:
            json.dump(output, f, indent=4) #added indent for readability
            
        # Write distances file
        with open(os.path.join(output_dir, "p_distances.txt"), "w") as f:
            for ref_id, distance in p_distances.items():
                f.write(f"{ref_id}: {distance}\n")
                
        # Save the combined alignment to a file
        output_alignment_path = os.path.join(output_dir, "updated_alignment.fasta")
        with open(output_alignment_path, "w") as f:
            SeqIO.write(alignment, output_alignment_path, "fasta")
        
        return output_alignment_path #returning path
    

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
