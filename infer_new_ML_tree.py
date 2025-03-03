import subprocess
import sys
import json

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    mafft_command = f"mafft --thread -1 --quiet --add {new_sequence} --keeplength {existing_alignment} > {output_alignment}"
    try:
        subprocess.run(mafft_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        return {"error": f"Failed to add sequence to alignment: {e}"}
    return None

def run_phylogenetic_placement(output_alignment, existing_tree):
    # Run IQ-TREE with the guide tree
    iqtree_command = f"iqtree2 -seed 2803 -nt 20 -redo --quiet -s {output_alignment} -g {existing_tree} -pre {output_alignment}_pp -m GTR+F+G4"
    try:
        subprocess.run(iqtree_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        return {"error": f"Failed to perform phylogenetic placement: {e}"}
    return None

def infer_global_optimization_tree(output_alignment, output_tree):
    # Run IQ-TREE with the constraint tree for optimization
    iqtree_command2 = f"iqtree2 -seed 2803 -nt 20 -redo --quiet -s {output_alignment} -t {output_alignment}_pp.treefile -pre {output_tree} -m GTR+F+G4"
    try:
        subprocess.run(iqtree_command2, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        return {"error": f"Failed to infer optimized tree: {e}"}
    return None

def main():
    if len(sys.argv) != 4 and len(sys.argv) != 6:
        print("Usage: python script_name.py existing_alignment.fasta new_sequence.fasta existing_tree.treefile [output_alignment.fasta output_tree_prefix]")
        sys.exit(1)
    
    existing_alignment = sys.argv[1]
    new_sequence = sys.argv[2]
    existing_tree = sys.argv[3]
    
    if len(sys.argv) == 6:
        output_alignment = sys.argv[4]
        output_tree = sys.argv[5]
    else:
        output_alignment = f"{existing_alignment.split('.')[0]}_updated.fasta"
        output_tree = f"{existing_tree.split('.')[0]}_reoptimised"
    
    # Add new sequence to existing alignment
    error = add_sequence_to_msa(existing_alignment, new_sequence, output_alignment)
    if error:
        with open("ml_tree_error.json", "w") as file:
            json.dump(error, file)
        sys.exit(1)
    
    # Run IQ-TREE phylogenetic placement first
    error = run_phylogenetic_placement(output_alignment, existing_tree)
    if error:
        with open("ml_tree_error.json", "w") as file:
            json.dump(error, file)
        sys.exit(1)
    
    # Run IQ-TREE with the constraint tree for optimization
    error = infer_global_optimization_tree(output_alignment, output_tree)
    if error:
        with open("ml_tree_error.json", "w") as file:
            json.dump(error, file)
        sys.exit(1)

    # Prepare output dictionary
    output = {
        "output_alignment": output_alignment,
        "output_tree": output_tree
    }

    # Write the output to a file
    with open("output/ml_tree_output.json", "w") as file:
        json.dump(output, file)

if __name__ == "__main__":
    main()
