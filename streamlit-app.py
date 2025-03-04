import streamlit as st
import subprocess
import sys
import os
import json
import tempfile
from Bio import SeqIO
import time
import zipfile

def log_error(message):
    st.error(message)

def extract_fasta_header(input_fasta):
    with open(input_fasta, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return record.id

def calculate_p_distance(input_fasta, reference_fasta):
    try:
        command = [sys.executable, "p-distance-calc.py", input_fasta, reference_fasta]
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"Error calculating p-distance. Return code: {result.returncode}")
            log_error(f"Standard Error:\n{result.stderr}")
            return None

        try:
            with open("output/p_distance_output.json", "r") as file:
                result = json.load(file)
            return result
        except FileNotFoundError:
            log_error("Error: p_distance_output.json not found in output directory.")
            return None
        except json.JSONDecodeError:
            log_error("Error: Could not decode JSON from p_distance_output.json")
            return None

    except subprocess.CalledProcessError as e:
        log_error(f"Subprocess error: {e}")
        return None

def infer_new_tree(existing_alignment, new_sequence, query_id, existing_tree, output_dir):
    output_alignment = os.path.join(output_dir, f"{query_id}_updated.fasta")
    output_tree = os.path.join(output_dir, f"{query_id}_reoptimised")
    script_path = "infer_new_ML_tree.py"

    command = [sys.executable, script_path, existing_alignment, new_sequence, existing_tree, output_alignment, output_tree]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"Error inferring new ML tree. Return code: {result.returncode}")
            log_error(f"Standard Error:\n{result.stderr}")
            return None, None, None

        try:
            with open("output/ml_tree_output.json", "r") as file:
                result = json.load(file)
            return result, output_alignment, output_tree
        except FileNotFoundError:
            log_error("Error: ml_tree_output.json not found in output directory.")
            return None, None, None
        except json.JSONDecodeError:
            log_error("Error: Could not decode JSON from ml_tree_output.json")
            return None, None, None

    except FileNotFoundError:
        log_error("Error: infer_new_ML_tree.py not found. Ensure the file exists and the path is correct.")
        return None, None, None
    except subprocess.CalledProcessError as e:
        log_error(f"Subprocess error: {e}")
        return None, None, None

def infer_subtype(input_newick, predefined_label, csv_file):
    command = [sys.executable, "ML_patristic-dist_calc.py", input_newick, predefined_label, csv_file]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"Error inferring subtype. Return code: {result.returncode}")
            log_error(f"Standard Error:\n{result.stderr}")
            return None

        try:
            with open("output/subtype_output.json", "r") as file:
                result = json.load(file)
            return result
        except FileNotFoundError:
            log_error("Error: subtype_output.json not found in output directory.")
            return None
        except json.JSONDecodeError:
            log_error("Error: Could not decode JSON from subtype_output.json")
            return None
    except subprocess.CalledProcessError as e:
        log_error(f"Subprocess error: {e}")
        return None

def main():
    st.title("Rat Hepatitis E Subtyping Tool v1.0")
    st.header("Sridhar Group")
    
    reference_fasta = "reference_genomes.fa"
    existing_alignment = "reference_alignment.fa"
    existing_tree = "reference_tree.tree"
    csv_file = "reference_subtypes.csv"

    input_fasta = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])

    output_dir = "output"
    
    if input_fasta is not None:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w+t") as tmp_file:
            tmp_file.write(input_fasta.getvalue().decode())
            temp_fasta_path = tmp_file.name

        query_id = extract_fasta_header(temp_fasta_path)

        st.write("\n**Summary Statistics:**")
        st.write(f"**Query ID:** {query_id}")
        
        with open(temp_fasta_path, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                st.write(f"**Query Length:** {len(record.seq)}")
                break
        
        os.makedirs(output_dir, exist_ok=True)
        
        progress_bar = st.progress(0)
        status_placeholder = st.empty()

        status_placeholder.write("\nCalculating p-distance...")
        for i in range(40):
            progress_bar.progress(i / 100)
            time.sleep(0.8)
        p_distance_output = calculate_p_distance(temp_fasta_path, reference_fasta)
        if p_distance_output:
            st.success("P-distance calculation completed.")
        else:
            st.error("Failed to calculate p-distance.")
            return
        
        status_placeholder.write("\nInferring new ML tree...")
        for i in range(59):
            progress_value = 0.4 + (i / 59) * 0.59
            progress_bar.progress(progress_value)
            time.sleep(0.8)
        tree_output, output_alignment, output_tree = infer_new_tree(existing_alignment, temp_fasta_path, query_id, existing_tree, output_dir)
        if tree_output:
            st.success("New ML tree inference completed.")
        else:
            st.error("Failed to infer new ML tree.")
            return
        
        status_placeholder.write("\nInferring subtype...")
        input_newick = f"{output_tree}.treefile"
        predefined_label = query_id
        subtype_output = infer_subtype(input_newick, predefined_label, csv_file)
        if subtype_output:
            st.success("Subtype inference completed.")
            progress_bar.progress(100)
        else:
            st.error("Failed to infer subtype.")
            return
        
        time.sleep(1)
        status_placeholder.write("\nAnalysis completed!")

        try:
            with open("output/p_distance_output.json", "r") as file:
                p_distance_result = json.load(file)
            with open("output/ml_tree_output.json", "r") as file:
                tree_result = json.load(file)
            with open("output/subtype_output.json", "r") as file:
                subtype_result = json.load(file)

            st.write("\n**Results Summary**")

            col1, col2 = st.columns(2)

            with col1:
                st.write("**P-Distance Results**")
                st.write(f"* **Closest Reference:** {p_distance_result['closest_reference']} ({p_distance_result['p_distance']:.4f})")
                st.write(f"* **Below Cutoff:** {p_distance_result['below_cutoff']}")

            with col2:
                st.write("**ML Patristic Distance Results**")
                st.write(f"* **Closest Reference:** {subtype_result['closest_reference_ml']} ({subtype_result['ml_distance']:.4f})")
                st.write(f"* **Conflicts:** {subtype_result['conflicts']}")

            if subtype_result['conflicts']:
                with st.expander("Conflict Summary"):
                    if subtype_result['subtype_assignment'] == "Not determined":
                        st.write("* **Consensus Assignment:** Not determined due to conflicts.")
                    else:
                        try:
                            conflicting_taxa_info = []
                            for taxon in subtype_result['conflict_summary']['conflicting_taxa']:
                                conflicting_taxa_info.append(f"{taxon['taxon']} (Clade {taxon['clade']} Subtype {taxon['subtype']})")
                            st.write(f"* **Conflicting Taxa:** {', '.join(conflicting_taxa_info)}")
                            st.write(f"* **Conflicting Clades:** {', '.join(subtype_result['conflict_summary']['clades'])}")
                            st.write(f"* **Conflicting Subtypes:** {', '.join(subtype_result['conflict_summary']['subtypes'])}")
                        except IndexError as e:
                            log_error(f"Error parsing subtype assignment: {e}")
            
            st.write(f"\n**Subtype Assignment:** {subtype_result['subtype_assignment']}")

        except FileNotFoundError as e:
            log_error(f"Error loading results: File not found - {e}")
        except json.JSONDecodeError as e:
            log_error(f"Error loading results: JSON decode error - {e}")
        except KeyError as e:
            log_error(f"Error loading results: Key error - {e}")
        except Exception as e:
            log_error(f"An unexpected error occurred: {e}")

        try:
            os.remove(temp_fasta_path)
        except FileNotFoundError:
            pass

        iqtree_dir = os.path.join(output_dir, "iqtree")
        os.makedirs(iqtree_dir, exist_ok=True)
        
        for filename in os.listdir(output_dir):
            if "reoptimised" in filename or "updated" in filename:
                file_path = os.path.join(output_dir, filename)
                destination_path = os.path.join(iqtree_dir, filename)
                os.rename(file_path, destination_path)
        
        with zipfile.ZipFile('output.zip', 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for root, dirs, files in os.walk('output'):
                for file in files:
                    file_path = os.path.join(root, file)
                    relative_path = os.path.relpath(file_path, start=os.path.dirname('output'))
                    zip_file.write(file_path, relative_path)

        with open("output.zip", "rb") as file:
            zip_data = file.read()

        st.download_button("Download Output", zip_data, file_name="output.zip")

if __name__ == "__main__":
    main()