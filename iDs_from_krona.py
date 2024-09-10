import re

# Path to your input and output files
input_file_path = 'kraken2_classified.txt'
output_file_path = 'seq_ids_and_taxids_from_krona.tsv'

# Define regex patterns to extract identifiers and taxids
pattern = re.compile(r'^>(\S+) .*kraken:taxid\|(\d+)', re.MULTILINE)

# Read the file and find matches
with open(input_file_path, 'r') as file:
    data = file.read()

# Find all matches
matches = pattern.findall(data)

# Write results to output file in TSV format
with open(output_file_path, 'w') as file:
    file.write('SeqID\tTaxID\n')  # Write header
    for seq_id, taxid in matches:
        file.write(f"{seq_id}\t{taxid}\n")

print(f"Results saved to {output_file_path}")
