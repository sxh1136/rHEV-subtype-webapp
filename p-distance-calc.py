import sys
import json
from Bio import SeqIO, Align

def calculate_p_distance(seq1, seq2):
    # Create a PairwiseAligner object with appropriate scoring
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"

    # Perform the alignment 
    alignments = aligner.align(seq1, seq2)

    alignment = next(alignments)

    # Extract aligned sequences as strings
    aligned_seq1 = str(alignment[0])
    aligned_seq2 = str(alignment[1])

    # Count mismatches and valid positions (ignoring gaps)
    mismatches = 0
    total_positions = 0
    
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != '-' or b != '-':  # Count all positions where at least one is valid
            total_positions += 1
            if a != b and a != '-' and b != '-':
                mismatches += 1
    
    # Calculate p-distance
    if total_positions == 0:
        return float('nan')  # Avoid division by zero if no valid positions
    p_distance = mismatches / total_positions
    return p_distance

def main(input_fasta, reference_fasta):
    try:
        input_seq = next(SeqIO.parse(input_fasta, "fasta")).seq
        input_id = next(SeqIO.parse(input_fasta, "fasta")).id  # Get the ID of the input sequence
        min_distance = float('inf')
        min_id = None
        
        # Read all reference sequences first to determine the total count
        reference_records = list(SeqIO.parse(reference_fasta, "fasta"))
        total_references = len(reference_records)
        
        # Dictionary to store p-distances for each reference
        p_distances = {}
        
        for index, record in enumerate(reference_records):
            ref_seq = record.seq
            distance = calculate_p_distance(str(input_seq), str(ref_seq))  
            
            # Store p-distance for this reference
            p_distances[record.id] = distance
            
            # Check if this distance is the smallest found so far
            if distance < min_distance:
                min_distance = distance
                min_id = record.id

        # Prepare output dictionary
        output = {
            "closest_reference": min_id,
            "p_distance": min_distance,
            "below_cutoff": min_distance <= 0.1833
        }

        # Write the output to a file
        with open("output/p_distance_output.json", "w") as file:
            json.dump(output, file)
        
        # Write p-distances to a file
        with open("output/p_distances.txt", "w") as file:
            for ref_id, distance in p_distances.items():
                file.write(f"{ref_id}: {distance}\n")

    except Exception as e:
        sys.stderr.write(f"An error occurred: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python p_distance.py <input_fasta> <reference_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    
    main(input_fasta, reference_fasta)
