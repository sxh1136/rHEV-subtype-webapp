import subprocess
import sys
import json
import os
import streamlit as st

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    mafft_command = f"mafft-linux64/mafft.bat --thread -1 --quiet --add {new_sequence} --keeplength {existing_alignment}"

    with open(output_alignment, 'w') as output_file:
        try:
            result = subprocess.run(mafft_command, check=True, stdout=output_file, stderr=subprocess.PIPE, text=True, shell=True)

            if result.returncode != 0:
                error_output = result.stderr.strip()
                st.error(f"MAFFT failed with return code {result.returncode}: {error_output}")
                return {"error": f"MAFFT failed with return code {result.returncode}: {error_output}"}

        except subprocess.CalledProcessError as e:
            st.error(f"MAFFT Error (CalledProcessError): {e}")
            return {"error": f"MAFFT Error (CalledProcessError): {e}"}
        except FileNotFoundError as e:
            st.error(f"MAFFT Error (FileNotFoundError): {e}")
            return {"error": f"MAFFT Error (FileNotFoundError): {e}"}

    if not os.path.exists(output_alignment):
        st.error(f"Output alignment file was not created: {output_alignment}")
        return {"error": f"Output alignment file was not created: {output_alignment}"}

    return None

def run_phylogenetic_placement(output_alignment, existing_tree):
    iqtree_command = f"./iqtree2 -seed 2803 -nt 20 -redo --quiet -s {output_alignment} -g {existing_tree} -pre {output_alignment}_pp -m GTR+F+G4"
    
    try:
        result = subprocess.run(iqtree_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)

        if result.returncode != 0:
            error_output = result.stderr.strip()
            st.error(f"Failed to perform phylogenetic placement: {error_output}")
            return {"error": f"Failed to perform phylogenetic placement: {error_output}"}
        return None
    except subprocess.CalledProcessError as e:
        st.error(f"Failed to perform phylogenetic placement: {e}")
        return {"error": f"Failed to perform phylogenetic placement: {e}"}

def infer_global_optimization_tree(output_alignment, output_tree):
    iqtree_command2 = f"./iqtree2 -seed 2803 -nt 20 -redo --quiet -s {output_alignment} -t {output_alignment}_pp.treefile -pre {output_tree} -m GTR+F+G4"
    
    try:
        result = subprocess.run(iqtree_command2, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)

        if result.returncode != 0:
            error_output = result.stderr.strip()
            st.error(f"Failed to infer optimized tree: {error_output}")
            return {"error": f"Failed to infer optimized tree: {error_output}"}
        return None
    except subprocess.CalledProcessError as e:
        st.error(f"Failed to infer optimized tree: {e}")
        return {"error": f"Failed to infer optimized tree: {e}"}

def main():
    if len(sys.argv) != 6:
        st.error("Usage: python script_name.py existing_alignment.fasta new_sequence.fasta existing_tree.treefile output_alignment.fasta output_tree_prefix")
        sys.exit(1)

    existing_alignment = os.path.abspath(sys.argv[1])
    new_sequence = os.path.abspath(sys.argv[2])
    existing_tree = os.path.abspath(sys.argv[3])
    output_alignment = os.path.abspath(sys.argv[4])
    output_tree = os.path.abspath(sys.argv[5])

    if not os.path.exists(existing_alignment):
        st.error(f"Existing alignment file not found: {existing_alignment}")
        sys.exit(1)
    if not os.path.exists(new_sequence):
        st.error(f"New sequence file not found: {new_sequence}")
        sys.exit(1)
    if not os.path.exists(existing_tree):
        st.error(f"Existing tree file not found: {existing_tree}")
        sys.exit(1)

    error = add_sequence_to_msa(existing_alignment, new_sequence, output_alignment)
    if error:
        with open("output/ml_tree_error.json", "w") as file:
            json.dump(error, file)
        st.error(f"Error during MAFFT: {error}")
        sys.exit(1)

    error = run_phylogenetic_placement(output_alignment, existing_tree)
    if error:
        with open("output/ml_tree_error.json", "w") as file:
            json.dump(error, file)
        st.error(f"Error during phylogenetic placement: {error}")
        sys.exit(1)

    error = infer_global_optimization_tree(output_alignment, output_tree)
    if error:
        with open("output/ml_tree_error.json", "w") as file:
            json.dump(error, file)
        st.error(f"Error during optimization tree inference: {error}")
        sys.exit(1)

    output = {
        "output_alignment": output_alignment,
        "output_tree": output_tree + ".treefile"
    }

    output_json_path = os.path.join(os.path.dirname(output_alignment), "ml_tree_output.json")
    with open(output_json_path, "w") as file:
        json.dump(output, file)

if __name__ == "__main__":
    main()