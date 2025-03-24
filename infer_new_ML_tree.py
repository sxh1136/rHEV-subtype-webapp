import subprocess
import sys
import json
import os
import streamlit as st
from Bio import SeqIO

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

    # Run IQ-TREE phylogenetic placement first
    error = run_phylogenetic_placement(existing_alignment, existing_tree) #no longer using the output alignment, but rather the existing_alignment which has already been aligned
    if error:
        error_file_path = os.path.join(output_dir, "ml_tree_error.json")
        with open(error_file_path, "w") as file:
            json.dump(error, file)
        print(f"Error during iqtree_command: {error}", file=sys.stderr)
        sys.exit(1)

    # Run IQ-TREE with the constraint tree for optimization
    error = infer_global_optimization_tree(existing_alignment, output_tree) #passing in the existing alignment file path
    if error:
        error_file_path = os.path.join(output_dir, "ml_tree_error.json")
        with open(error_file_path, "w") as file:
            json.dump(error, file)
        print(f"Error during iqtree_command2: {error}", file=sys.stderr)
        sys.exit(1)

    # Prepare output dictionary
    output = {
        "output_alignment": existing_alignment,
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
