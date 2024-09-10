# Define file paths
seqid_file = 'seqID_species_taxoID.txt'
ph_file = 'extracted_seqIDs_from_ph_file.txt'
output_file = 'dictionary_for_tree.txt'  # Optionally, you can overwrite the original file

# Read the seqID_species_taxoID.txt file and create a dictionary
seqid_dict = {}
with open(seqid_file, 'r') as file:
    header = file.readline()  # Skip header if present
    for line in file:
        parts = line.strip().split('\t')  # Adjust delimiter if needed
        seq_id, taxon_id, species_name = parts
        seqid_dict[seq_id] = (taxon_id, species_name)

# Read the extracted_seqIDs_from_ph_file.txt file and add Taxon_ID and Species_Name
with open(ph_file, 'r') as file:
    lines = file.readlines()

# Assume the first line is the header; otherwise, adjust accordingly
header = lines[0].strip().split('\t')
header.extend(['Taxon_ID', 'Species_Name'])  # Add new column titles
updated_lines = [header]

for line in lines[1:]:
    parts = line.strip().split('\t')  # Adjust delimiter if needed
    seq_id = parts[0]  # Assuming seqID is the first column; adjust if necessary
    if seq_id in seqid_dict:
        taxon_id, species_name = seqid_dict[seq_id]
        parts.extend([taxon_id, species_name])
    else:
        parts.extend(['', ''])  # Add empty fields if no match
    updated_lines.append(parts)

# Write the updated data to the new file
with open(output_file, 'w') as file:
    for line in updated_lines:
        file.write('\t'.join(line) + '\n')

print(f'Updated file saved as {output_file}')