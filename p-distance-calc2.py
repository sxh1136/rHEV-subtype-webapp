import sys
import json
import os
import tempfile
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import subprocess
import numpy as np
os.environ['R_HOME'] = '/usr/lib/R'
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import pandas as pd

# Activate pandas conversion
pandas2ri.activate()


def create_distance_matrix(alignment):
    """
    Calculates a distance matrix using the K80 model with Gamma correction (gamma=4)
    using R's ape package.

    Args:
        alignment (MultipleSeqAlignment): A BioPython MultipleSeqAlignment object.

    Returns:
        dict: A dictionary representing the distance matrix, where keys are tuples of sequence IDs
              and values are the K80 distances.
    """
    ids = [seq.id for seq in alignment]

    # Convert alignment to FASTA string
    fasta_string = ""
    for seq in alignment:
        fasta_string += f">{seq.id}\n{str(seq.seq)}\n"

    # Use R to calculate the distance matrix
    distance_matrix, base_frequencies = calculate_k80_distance_rpy2(fasta_string)

    # Convert the R matrix to a Python dictionary (though it's not really used as such)
    matrix = {}
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            matrix[(ids[i], ids[j])] = distance_matrix[i, j]
            matrix[(ids[j], ids[i])] = distance_matrix[j, i]  # Matrix is symmetric
    return matrix, base_frequencies


def calculate_k80_distance_rpy2(fasta_string):
    """
    Calculates the K80 genetic distance matrix using R's ape package via rpy2.

    Args:
        fasta_string (str): A string containing the alignment in FASTA format.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: A NumPy array representing the K80 distance matrix.
            - dict: A dictionary of base frequencies.
    """
    try:
        # Import R packages
        ape = importr('ape')
        msa = importr('msa')

        # Define R code as a string
        r_code = """
        require(ape)
        require(msa)
        require(stringr)

        # Function to calculate K80 distance matrix and base frequencies
        calculate_distance_matrix <- function(fasta_string) {
            try {
                # Create a temporary file to store the FASTA string
                temp_file <- tempfile(fileext = ".fasta")
                writeLines(fasta_string, temp_file)

                # Read the alignment from the temporary file
                aln <- read.alignment(file = temp_file, format = "fasta")

                # Check if the alignment is valid
                if (is.null(aln$seq) || length(aln$seq) == 0) {
                    stop("Alignment file is empty or invalid.")
                }

                # Convert the alignment to a DNAbin object
                aln.2 <- msaConvert(aln, 'ape::DNAbin')

                # Calculate the K80 distance matrix
                dist_matrix <- dist.dna(aln.2, model = "K80", gamma = 4, as.matrix = TRUE)

                # Clean up the temporary file
                unlink(temp_file)

                # Estimate base frequencies
                frequencies <- base.freq(aln.2)

                result <- list(dist_matrix = dist_matrix, frequencies = frequencies)

                return(result)

            } catch (e) {
                # Handle errors within R
                print(paste("R Error:", e$message))
                return(NULL)  # Return NULL to indicate failure
            }
        }

        # Call the function and return the result
        result <- calculate_distance_matrix(fasta_string=fasta_string)
        """


        # Execute the R code
        robjects.r(r_code)

        # Extract the result from R
        r_result = robjects.r['result']

        if r_result is robjects.NULL:
            raise ValueError("R function returned NULL, indicating an error.")

        distance_matrix = np.array(r_result[0])
        frequencies = dict(zip(robjects.r('names(result$frequencies)'), r_result[1])) #Extract base frequencies

        return distance_matrix, frequencies

    except Exception as e:
        print(f"Python Error: {e}")
        raise  # Re-raise the exception to be caught in the main function


def main(input_fasta, reference_fasta, output_dir):
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
            muscle_cline = ["./muscle", "-profile", "-in1", reference_fasta, "-in2", new_seq_file, "-out", combined_alignment_file]
            subprocess.run(muscle_cline, check=True, capture_output=True, text=True)
            
            # Read the combined alignment
            aligned_sequences = list(SeqIO.parse(combined_alignment_file, "fasta"))
        
        # Create MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(aligned_sequences)

        
        # Calculate distance matrix using K80 model with gamma correction and base frequencies
        distance_matrix, base_frequencies = create_distance_matrix(alignment)
        
        # Extract distances for input sequence
        p_distances = {}
        ids = [seq.id for seq in alignment]

        for i, seq in enumerate(aligned_sequences):
            if seq.id != new_input_id:
                ref_id = seq.id
                # Access distance from the matrix, order matters!
                ref_index = ids.index(ref_id)
                query_index = ids.index(new_input_id)
                distance = distance_matrix[query_index, ref_index]

                p_distances[ref_id] = distance
        
        # Find closest reference
        min_id, min_distance = min(p_distances.items(), key=lambda x: x[1])
        
        # Prepare output
        output = {
            "closest_reference": min_id,
            "p_distance": min_distance,
            "below_cutoff": min_distance <= 0.1833,
            "original_input_id": original_id,  # Include the original ID in the output
            "gamma": 4.0, #Hardcoded gamma value
            "base_frequencies": base_frequencies  # Include base frequencies in the output
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
        print("Usage: python msa_pdistance.py <input_fasta> <reference_fasta> <output_dir>")
        sys.exit(1)
        
    alignment_path = main(sys.argv[1], sys.argv[2], sys.argv[3])