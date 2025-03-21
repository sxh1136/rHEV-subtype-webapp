import subprocess
import sys
import json
import os
import streamlit as st

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    # Construct the MAFFT command as a string
    mafft_command = f"mafft-linux64/mafft.bat --thread -1 --quiet --add {new_sequence} --keeplength {existing_alignment}"

    with open(output_alignment, 'w') as output_file:
        try:
            # Run the MAFFT command with shell=True
            result = subprocess.run(mafft_command, check=True, stdout=output_file, stderr=subprocess.PIPE, text=True, shell=True)

            if result.returncode != 0:
                error_output = result.stderr.strip()
                print(f"MAFFT Error: MAFFT failed with return code {result.returncode}", file=sys.stderr)
                print(f"MAFFT Error (stderr):\n{error_output}", file=sys.stderr)
                return {"error": f"MAFFT failed with return code {result.returncode}: {error_output}"}

        except subprocess.CalledProcessError as e:
            print(f"MAFFT Error (CalledProcessError): {e}", file=sys.stderr)
            return {"error": f"MAFFT Error (CalledProcessError): {e}"}
        except FileNotFoundError as e:
            print(f"MAFFT Error (FileNotFoundError): {e}", file=sys.stderr)
            return {"error": f"MAFFT Error (FileNotFoundError): {e}"}

    # Check that the output_alignment was created:
    if not os.path.exists(output_alignment):
        print(f"ERROR: add_sequence_to_msa: Output alignment file was not created: {output_alignment}", file=sys.stderr)
        return {"error": f"add_sequence_to_msa: Output alignment file was not created: {output_alignment}"}

    return None

def log_file_contents(file_path):
    if os.path.exists(file_path):
        with open(file_path) as f:
            st.write(f"Contents of {file_path}:\n", f.readlines()[:5])  # Log the first few lines
    else:
        st.error(f"File not found: {file_path}")


def run_phylogenetic_placement(output_alignment, existing_tree):
    # Construct the IQ-TREE command as a string
    iqtree_command = f"./iqtree2 -seed 2803 -redo -s {output_alignment} -g {existing_tree} -pre {output_alignment}_pp -m GTR+F+G4"
    
    st.write(f"Running command: {iqtree_command}")  # Log the command

    try:
        result = subprocess.run(iqtree_command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        st.write("IQ-TREE Standard Output:", result.stdout)
        st.write("IQ-TREE Standard Error:", result.stderr)

        if result.returncode != 0:
            error_output = result.stderr.strip()
            st.error(f"Failed to perform phylogenetic placement: {error_output}")
            return {"error": f"Failed to perform phylogenetic placement: {error_output}"}
    except Exception as e:
        st.error(f"Exception during IQ-TREE execution: {str(e)}")
        return {"error": f"Exception during IQ-TREE execution: {str(e)}"}

    return None

def infer_global_optimization_tree(output_alignment, output_tree):
    # Construct the IQ-TREE command for optimization as a string
    iqtree_command2 = f"./iqtree2 -seed 2803 -redo -s {output_alignment} -t {output_alignment}_pp.treefile -pre {output_tree} -m GTR+F+G4"

    try:
        result = subprocess.run(iqtree_command2, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        st.write(result.stdout)  # Log standard output
        st.write(result.stderr)   # Log standard error
        if result.returncode != 0:
            error_output = result.stderr.strip()
            print(f"iqtree_command2 Error: {error_output}", file=sys.stderr)
            return {"error": f"Failed to infer optimized tree: {error_output}"}
        return None
    except subprocess.CalledProcessError as e:
        print(f"iqtree_command2 Error: {e}", file=sys.stderr)
        return {"error": f"Failed to infer optimized tree: {e}"}

def main(existing_alignment, new_sequence, existing_tree, output_alignment, output_tree, output_dir):
    # Check file existence *before* doing anything else
    if not os.path.exists(existing_alignment):
        print(f"ERROR: Existing alignment file not found: {existing_alignment}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(new_sequence):
        print(f"ERROR: New sequence file not found: {new_sequence}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(existing_tree):
        print(f"ERROR: Existing tree file not found: {existing_tree}", file=sys.stderr)
        sys.exit(1)

    # Add new sequence to existing alignment
    error = add_sequence_to_msa(existing_alignment, new_sequence, output_alignment)
    if error:
        error_file_path = os.path.join(output_dir, "ml_tree_error.json")
        with open(error_file_path, "w") as file:
            json.dump(error, file)
        print(f"Error during mafft: {error}", file=sys.stderr)
        sys.exit(1)

    # Call this function for the necessary files
    log_file_contents(existing_alignment)
    log_file_contents(output_alignment)  # Check the output alignment file
    log_file_contents(existing_tree)

    # Run IQ-TREE phylogenetic placement first
    error = run_phylogenetic_placement(output_alignment, existing_tree)
    if error:
        error_file_path = os.path.join(output_dir, "ml_tree_error.json")
        with open(error_file_path, "w") as file:
            json.dump(error, file)
        print(f"Error during iqtree_command: {error}", file=sys.stderr)
        sys.exit(1)

    # Run IQ-TREE with the constraint tree for optimization
    error = infer_global_optimization_tree(output_alignment, output_tree)
    if error:
        error_file_path = os.path.join(output_dir, "ml_tree_error.json")
        with open(error_file_path, "w") as file:
            json.dump(error, file)
        print(f"Error during iqtree_command2: {error}", file=sys.stderr)
        sys.exit(1)

    # Prepare output dictionary
    output = {
        "output_alignment": output_alignment,
        "output_tree": output_tree + ".treefile"
    }

    # Write the output to a file
    output_json_path = os.path.join(output_dir, "ml_tree_output.json")
    with open(output_json_path, "w") as file:
        json.dump(output, file)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script_name.py existing_alignment.fasta new_sequence.fasta existing_tree.treefile output_alignment.fasta output_tree_prefix output_dir", file=sys.stderr)
        sys.exit(1)

    existing_alignment = os.path.abspath(sys.argv[1])
    new_sequence = os.path.abspath(sys.argv[2])
    existing_tree = os.path.abspath(sys.argv[3])
    output_alignment = os.path.abspath(sys.argv[4])
    output_tree = os.path.abspath(sys.argv[5])
    output_dir = os.path.abspath(sys.argv[6])

    main(existing_alignment, new_sequence, existing_tree, output_alignment, output_tree, output_dir)
