import re

# Define file paths
input_file = 'kraken2_classified.txt'
output_file = 'krona_classified.fasta'

# Open the input file and the output FASTA file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):  # Line with sequence ID and metadata
            # Write the entire header line in FASTA format
            outfile.write(line)  # Write the header line as is
        else:  # Line with sequence data
            # Write the sequence data
            outfile.write(line)

print(f'FASTA file created: {output_file}')
