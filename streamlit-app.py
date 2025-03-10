import streamlit as st
import subprocess
import sys
import os
import json
import tempfile
from Bio import SeqIO
import time
import zipfile
from datetime import datetime
import pandas as pd

def log_error(message):
    st.error(message)

def extract_fasta_header(input_fasta):
    try:
        with open(input_fasta, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                return record.id
    except Exception as e:
        log_error(f"Error extracting FASTA header: {e}")
        return None

def calculate_p_distance(input_fasta, reference_fasta, output_dir):
    try:
        command = [sys.executable, "p-distance-calc.py", input_fasta, reference_fasta, output_dir]
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"Error calculating p-distance. Return code: {result.returncode}")
            log_error(f"Standard Error:\n{result.stderr}")
            return None

        try:
            with open(os.path.join(output_dir, "p_distance_output.json"), "r") as file:
                result = json.load(file)
            return result
        except FileNotFoundError:
            log_error(f"Error: p_distance_output.json not found in {output_dir} directory.")
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

    command = [sys.executable, script_path, existing_alignment, new_sequence, existing_tree, output_alignment, output_tree, output_dir]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            st.error(f"Error inferring new ML tree. Return code: {result.returncode}")
            return None, None, None

        # Check if the output file was created
        ml_tree_output_path = os.path.join(output_dir, "ml_tree_output.json")
        if not os.path.exists(ml_tree_output_path):
            st.error(f"Error: {ml_tree_output_path} not found after execution.")
            return None, None, None

        with open(ml_tree_output_path, "r") as file:
            result = json.load(file)
        return result, output_alignment, output_tree

    except FileNotFoundError:
        st.error("Error: infer_new_ML_tree.py not found. Ensure the file exists and the path is correct.")
        return None, None, None
    except subprocess.CalledProcessError as e:
        st.error(f"Subprocess error: {e}")
        st.error(f"Return code: {e.returncode}")
        st.error(f"Output: {e.output}")
        st.error(f"Error: {e.stderr}")
        return None, None, None

def infer_subtype(input_newick, predefined_label, csv_file, output_dir):
    command = [sys.executable, "ML_patristic-dist_calc.py", input_newick, predefined_label, csv_file, output_dir]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"Error inferring subtype. Return code: {result.returncode}")
            log_error(f"Standard Error:\n{result.stderr}")
            return None

        try:
            with open(os.path.join(output_dir, "subtype_output.json"), "r") as file:
                result = json.load(file)
            return result
        except FileNotFoundError:
            log_error(f"Error: subtype_output.json not found in {output_dir} directory.")
            return None
        except json.JSONDecodeError:
            log_error("Error: Could not decode JSON from subtype_output.json")
            return None
    except subprocess.CalledProcessError as e:
        log_error(f"Subprocess error: {e}")
        return None

def main():
    st.image("header.jpg")
    st.title("Rat Hepatitis E Subtyping Tool v1.0")
    st.header("Sridhar Group")
    
    reference_fasta = "reference_genomes.fa"
    existing_alignment = "reference_alignment.fa"
    existing_tree = "reference_tree.tree"
    csv_file = "reference_subtypes.csv"

    # Input options: file upload or text area
    input_option = st.radio("Choose input method:", ("Upload FASTA file", "Paste FASTA sequence"))
    if input_option == "Upload FASTA file":
        input_fasta = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])
        fasta_string = None  # Initialize to None
    else:
        fasta_string = st.text_area("Paste FASTA sequence", height=200)
        input_fasta = None  # Initialize to None

    # Add a button to start the analysis
    start_analysis = st.button("Start Analysis")

    if start_analysis:
        if input_fasta:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w+t") as tmp_file:
                tmp_file.write(input_fasta.getvalue().decode())
                temp_fasta_path = tmp_file.name
        elif fasta_string: #Using the text area input
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w+t") as tmp_file:
                tmp_file.write(fasta_string)
                temp_fasta_path = tmp_file.name
        else:
            st.warning("Please upload a FASTA file or paste a FASTA sequence.")
            return

        query_id = extract_fasta_header(temp_fasta_path)
        if query_id is None:
            st.error("Could not extract FASTA header. Please check your input.")
            try:
                os.remove(temp_fasta_path)
            except FileNotFoundError:
                pass
            return

        st.write("\n**Summary Statistics:**")
        st.write(f"**Query ID:** {query_id}")
        
        try:
            with open(temp_fasta_path, 'r') as file:
                for record in SeqIO.parse(file, "fasta"):
                    st.write(f"**Query Length:** {len(record.seq)}")
                    break
        except Exception as e:
            log_error(f"Error reading FASTA file: {e}")
            try:
                os.remove(temp_fasta_path)
            except FileNotFoundError:
                pass
            return
        
        # Generate dynamic output directory name
        now = datetime.now()
        timestamp = now.strftime("%Y%m%d")
        output_dir = f"rHEV_subtyping_results_{query_id}_{timestamp}"
        os.makedirs(output_dir, exist_ok=True)
        
        progress_bar = st.progress(0)
        status_placeholder = st.empty()
        status_placeholder.write("\nCalculating p-distance...")
        for i in range(40):
            progress_bar.progress(i / 100)
            time.sleep(0.3)
        p_distance_output = calculate_p_distance(temp_fasta_path, reference_fasta, output_dir)
        if p_distance_output:
            st.success("P-distance calculation completed.")
        else:
            st.error("Failed to calculate p-distance.")
            return
        
        status_placeholder.write("\nInferring new ML tree...")
        for i in range(59):
            progress_value = 0.4 + (i / 59) * 0.59
            progress_bar.progress(progress_value)
            time.sleep(0.3)
        tree_output, output_alignment, output_tree = infer_new_tree(existing_alignment, temp_fasta_path, query_id, existing_tree, output_dir)
        if tree_output:
            st.success("New ML tree inference completed.")
        else:
            st.error("Failed to infer new ML tree.")
            return
        
        status_placeholder.write("\nInferring subtype...")
        input_newick = f"{output_tree}.treefile"
        predefined_label = query_id
        subtype_output = infer_subtype(input_newick, predefined_label, csv_file, output_dir)
        if subtype_output:
            st.success("Subtype inference completed.")
            progress_bar.progress(100)
        else:
            st.error("Failed to infer subtype.")
            return
        
        time.sleep(1)
        status_placeholder.write("\nAnalysis completed!")

        try:
            with open(os.path.join(output_dir, "p_distance_output.json"), "r") as file:
                p_distance_result = json.load(file)
            with open(os.path.join(output_dir, "ml_tree_output.json"), "r") as file:
                tree_result = json.load(file)
            with open(os.path.join(output_dir, "subtype_output.json"), "r") as file:
                subtype_result = json.load(file)

            st.write("\n**Results Summary**")

            # Prepare data for the table, making metrics the column titles
            table_data = {
                "Closest Reference by P-Distance (value)": [
                    f"{p_distance_result['closest_reference']} ({p_distance_result['p_distance']:.4f})"
                ],
                "Closest Reference by ML Patristic Distance (value)": [
                    f"{subtype_result['closest_reference_ml']} ({subtype_result['ml_distance']:.4f})"
                ],
                "Clade Assignment": [
                    subtype_result['clade_assignment']
                ],
                "Subtype Assignment": [
                    subtype_result['subtype_assignment']
                ]
            }
        
            # Create a Pandas DataFrame
            df = pd.DataFrame(table_data)

            # Display the table using Streamlit
            st.table(df)
            st.write("**Note**: clade and subtype assignments are only determined if ML patristic distance is below cutoff values.")
            
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
        
        zip_filename = f"{output_dir}.zip"
        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    relative_path = os.path.relpath(file_path, start=output_dir)
                    zip_file.write(file_path, relative_path)

        with open(zip_filename, "rb") as file:
            zip_data = file.read()

        st.download_button("Download Output", zip_data, file_name=zip_filename)

        st.balloons()
if __name__ == "__main__":
    main()
