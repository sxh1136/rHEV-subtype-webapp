import sys
import dendropy
import os
import csv
import json

def main(input_newick, predefined_label, csv_file, output_dir):
    try:
        # Check if the input files exist
        if not os.path.isfile(input_newick):
            print(f"Error: File '{input_newick}' does not exist.")
            return
        if not os.path.isfile(csv_file):
            print(f"Error: File '{csv_file}' does not exist.")
            return

        # Load the tree from the Newick file
        tree = dendropy.Tree.get(
            path=input_newick,
            schema="newick",
            preserve_underscores=True
        )
        print("Tree loaded successfully.")

        # Create a phylogenetic distance matrix
        pdc = tree.phylogenetic_distance_matrix()
        print("Distance matrix created.")

         # Print the distance matrix
        print("Distance Matrix:")
        for taxon1 in tree.taxon_namespace:
            for taxon2 in tree.taxon_namespace:
                if taxon1 != taxon2:
                    distance = pdc(taxon1, taxon2)
                    print(f"Distance from {taxon1.label} to {taxon2.label}: {distance}")


        # Find the predefined taxon
        predefined_taxon = None
        for taxon in tree.taxon_namespace:
            if taxon.label == predefined_label:
                predefined_taxon = taxon
                break

        if predefined_taxon is None:
            print(f"Error: Taxon '{predefined_label}' not found in the tree.")
            return

        # List to store distances and corresponding taxa
        distance_list = []

        # Calculate distances from the predefined taxon to all other taxa
        for taxon in tree.taxon_namespace:
            if taxon != predefined_taxon:  # Skip the predefined taxon itself
                distance = pdc(predefined_taxon, taxon)
                print(f"Distance from {predefined_taxon.label} to {taxon.label}: {distance}")
                distance_list.append((distance, taxon.label))

        print(f"Distance list populated with {len(distance_list)} entries.")

        # Sort the list by distance (first element of tuple)
        distance_list.sort(key=lambda x: x[0])
        print("Distance list sorted.")

        print("Sorted distance list:")
        for dist, label in distance_list:
            print(f"{label}: {dist}")

        # ML nucleotide threshold
        clade_threshold = 1.0227
        subtype_threshold = 0.5586

        # Get the closest reference
        closest_distance, closest_label = distance_list[0]

        # Read the CSV file
        csv_data = {}
        with open(csv_file, mode='r', encoding='utf-8-sig') as file:
            reader = csv.DictReader(file, delimiter=',')
            for row in reader:
                csv_data[row['Name']] = {'Clade': row['Clade'], 'Subtype': row['Subtype']}

        print("CSV data loaded.")

        # Determine clade and subtype of the closest reference
        clade_assignment = "Not-determined"
        subtype_assignment = "Not-determined"
        if closest_label in csv_data:
            if closest_distance <= clade_threshold:
                clade_assignment = csv_data[closest_label]['Clade']
            if closest_distance <= subtype_threshold:
                subtype_assignment = csv_data[closest_label]['Subtype']

        print("Subtype and clade assignments determined.")

        # Prepare output dictionary
        output = {
            "closest_reference_ml": closest_label,
            "ml_distance": closest_distance,
            "below_cutoff": closest_distance < subtype_threshold,
            "clade_assignment": clade_assignment,
            "subtype_assignment": subtype_assignment
        }

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Write the output to a file in the specified directory
        output_file = os.path.join(output_dir, "subtype_output.json")
        with open(output_file, "w") as file:
            json.dump(output, file)

        print(f"Output written to {output_file}.")

        # Write patristic distances to a file
        distances_file = os.path.join(output_dir, "patristic_distances.txt")
        with open(distances_file, "w") as file:
            file.write(f"Patristic Distances from {predefined_label}:\n")
            for dist, label in distance_list:
                file.write(f"{label}: {dist}\n")

        print(f"Patristic distances written to {distances_file}.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_newick> <predefined_label> <csv_file> <output_dir>")
        sys.exit(1)

    input_newick = sys.argv[1]
    predefined_label = sys.argv[2]
    csv_file = sys.argv[3]
    output_dir = sys.argv[4]

    main(input_newick, predefined_label, csv_file, output_dir)
