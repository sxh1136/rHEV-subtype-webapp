import subprocess
import sys
import json
import os
import streamlit as st

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    mafft_command = ["./mafft-linux64/mafft.bat", "--thread", "-1", "--quiet", "--add", new_sequence, "--keeplength", existing_alignment]
    
    with open(output_alignment, 'w') as output_file:
        try:
            result = subprocess.run(mafft_command, check=True, stdout=output_file, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                error_output = result.stderr.strip()
                st.error(f"MAFFT Command Error: {error_output}")
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
    iqtree_command = [
        "./iqtree2", "-seed", "2803", "-nt", "20", "-redo", 
        "--quiet", "-s", output_alignment, "-g", existing_tree, 
        "-pre", f"{output_alignment}_pp", "-m", "GTR+F+G4"
    ]
    
    try:
        result = subprocess.run(iqtree_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            error_output = result.stderr.strip()
            st.error(f"IQ-TREE Command: {' '.join(iqtree_command)}")
            st.error(f"IQ-TREE Error Output: {error_output}")
            return {"error": f"Failed to perform phylogenetic placement: {error_output}"}
        return None
    except subprocess.CalledProcessError as e:
        st.error(f"IQ-TREE Error (CalledProcessError): {e}")
        return {"error": f"Failed to perform phylogenetic placement: {e}"}
    except Exception as e:
        st.error(f"IQ-TREE Unexpected Error: {e}")
        return {"error": f"An unexpected error occurred: {e}"}

def infer_global_optimization_tree(output_alignment, output_tree):
    iqtree_command2 = [
        "./iqtree2", "-seed", "2803", "-nt", "20", "-redo", 
        "--quiet", "-s", output_alignment, "-t", f"{output_alignment}_pp.treefile", 
        "-pre", output_tree, "-m", "GTR+F+G4"
    ]
    
    try:
        result = subprocess.run(iqtree_command2, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            error_output = result.stderr.strip()
            st.error(f"IQ-TREE Optimization Command: {' '.join(iqtree_command2)}")
            st.error(f"IQ-TREE Optimization Error Output: {error_output}")
            return {"error": f"Failed to infer optimized tree: {error_output}"}
        return None
    except subprocess.CalledProcessError as e:
        st.error(f"IQ-TREE Optimization Error (CalledProcessError): {e}")
        return {"error": f"Failed to infer optimized tree: {e}"}
    except Exception as e:
        st.error(f"IQ-TREE Optimization Unexpected Error: {e}")
        return {"error": f"An unexpected error occurred: {e}"}

def main():
    st.title("Phylogenetic Analysis Tool")

    if st.button("Run Analysis"):
        existing_alignment = st.text_input("Existing Alignment File")
        new_sequence = st.text_input("New Sequence File")
        existing_tree = st.text_input("Existing Tree File")
        output_alignment = st.text_input("Output Alignment File")
        output_tree = st.text_input("Output Tree Prefix")

        # Ensure all inputs are provided
        if not all([existing_alignment, new_sequence, existing_tree, output_alignment, output_tree]):
            st.error("Please provide all input files.")
            return

        # Check file existence
        if not os.path.exists(existing_alignment):
            st.error(f"Existing alignment file not found: {existing_alignment}")
            return
        if not os.path.exists(new_sequence):
            st.error(f"New sequence file not found: {new_sequence}")
            return
        if not os.path.exists(existing_tree):
            st.error(f"Existing tree file not found: {existing_tree}")
            return

        # Add new sequence to existing alignment
        error = add_sequence_to_msa(existing_alignment, new_sequence, output_alignment)
        if error:
            st.json(error)
            return

        # Run IQ-TREE phylogenetic placement
        error = run_phylogenetic_placement(output_alignment, existing_tree)
        if error:
            st.json(error)
            return

        # Run IQ-TREE with the constraint tree for optimization
        error = infer_global_optimization_tree(output_alignment, output_tree)
        if error:
            st.json(error)
            return

        # Prepare output dictionary
        output = {
            "output_alignment": output_alignment,
            "output_tree": output_tree + ".treefile"
        }

        # Write the output to a file
        output_json_path = os.path.join(os.path.dirname(output_alignment), "ml_tree_output.json")
        with open(output_json_path, "w") as file:
            json.dump(output, file)

        st.success("Analysis completed successfully!")

if __name__ == "__main__":
    main()