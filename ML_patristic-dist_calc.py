import sys
import dendropy
import os
import csv
from collections import Counter
import json

def main(input_newick, predefined_label, csv_file):
    # Check if the input files exist
    if not os.path.isfile(input_newick):
        return {"error": f"File '{input_newick}' does not exist."}
    if not os.path.isfile(csv_file):
        return {"error": f"File '{csv_file}' does not exist."}

    try:
        # Load the tree from the Newick file
        tree = dendropy.Tree.get(
            path=input_newick,
            schema="newick"
        )
    except Exception as e:
        return {"error": f"Failed to load tree from '{input_newick}': {e}"}

    # Create a phylogenetic distance matrix
    pdc = tree.phylogenetic_distance_matrix()

    # Find the predefined taxon
    predefined_taxon = None
    for taxon in tree.taxon_namespace:
        if taxon.label == predefined_label:
            predefined_taxon = taxon
            break

    if predefined_taxon is None:
        return {"error": f"Taxon '{predefined_label}' not found in the tree."}

    # List to store distances and corresponding taxa
    distance_list = []

    # Calculate distances from the predefined taxon to all other taxa
    for taxon in tree.taxon_namespace:
        if taxon != predefined_taxon:  # Skip the predefined taxon itself
            distance = pdc(predefined_taxon, taxon)
            distance_list.append((distance, taxon.label))

    # Sort the list by distance (first element of tuple)
    distance_list.sort(key=lambda x: x[0])

    # ML nucleotide threshold
    threshold = 0.5585

    # Filter distances under the threshold
    filtered_distances = [(dist, label) for dist, label in distance_list if dist < threshold]

    # Read the CSV file
    csv_data = {}
    with open(csv_file, mode='r', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file, delimiter=',')
        for row in reader:
            csv_data[row['Name']] = {'Clade': row['Clade'], 'Subtype': row['Subtype']}

    # Check for conflicting clades and subtypes
    clades = []
    subtypes = []
    conflicting_taxa_info = []
    for dist, label in filtered_distances:
        if label in csv_data:
            clades.append(csv_data[label]['Clade'])
            subtypes.append(csv_data[label]['Subtype'])
            if len(set(clades)) > 1 or len(set(subtypes)) > 1:
                conflicting_taxa_info.append({
                    'taxon': label,
                    'clade': csv_data[label]['Clade'],
                    'subtype': csv_data[label]['Subtype']
                })

    conflicts = False
    most_common_clade = None
    most_common_subtype = None
    if clades and subtypes:
        if len(set(clades)) > 1 or len(set(subtypes)) > 1:
            conflicts = True
            # In case of conflicts, you could still provide a most common clade/subtype
            # or handle it as "Not determined" in your output
            most_common_clade = Counter(clades).most_common(1)[0][0]
            most_common_subtype = Counter(subtypes).most_common(1)[0][0]
        else:
            most_common_clade = Counter(clades).most_common(1)[0][0]
            most_common_subtype = Counter(subtypes).most_common(1)[0][0]

    # Extract closest label after processing filtered distances
    closest_label = None
    if filtered_distances:  # Only try to get the label if filtered_distances is not empty
        closest_label = filtered_distances[0][1]  # First element's label

    # Prepare output dictionary
    output = {
        "closest_reference_ml": closest_label,
        "ml_distance": filtered_distances[0][0] if filtered_distances else None,
        "conflicts": conflicts,
        "conflict_summary": {
            "conflicting_taxa": [{"taxon": taxon['taxon'], "clade": taxon['clade'], "subtype": taxon['subtype']} for taxon in conflicting_taxa_info],
            "clades": list(set(clades)),
            "subtypes": list(set(subtypes))
        },
        "subtype_assignment": f"Clade={most_common_clade}, Subtype={most_common_subtype}" if most_common_clade and most_common_subtype else "Clade=Unknown, Subtype=Unknown"
    }

    # Write the output to a file
    with open("output/subtype_output.json", "w") as file:
        json.dump(output, file)

    # Write patristic distances to a file
    with open("output/patristic_distances.txt", "w") as file:
        file.write(f"Patristic Distances from {predefined_label}:\n")
        for dist, label in distance_list:
            file.write(f"{label}: {dist}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_newick> <predefined_label> <csv_file>")
        sys.exit(1)

    input_newick = sys.argv[1]
    predefined_label = sys.argv[2]
    csv_file = sys.argv[3]
    
    main(input_newick, predefined_label, csv_file)
