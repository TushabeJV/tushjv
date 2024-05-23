import subprocess
import os

# Define the input files
#R1_files = ["Trimmed_paired/BO4_trimmed_r1_paired.fq.gz", "Trimmed_paired/B05_trimmed_r1_paired.fq.gz", "Trimmed_paired/G03_trimmed_r1_paired.fq.gz", "Trimmed_paired/G11_trimmed_r1_paired.fq.gz"]
#R2_files = ["Trimmed_paired/BO4_trimmed_r2_paired.fq.gz", "Trimmed_paired/B05_trimmed_r2_paired.fq.gz", "Trimmed_paired/G03_trimmed_r2_paired.fq.gz", "Trimmed_paired/G11_trimmed_r2_paired.fq.gz"]

# Set the output directory to a new folder called "assembly_output" in the current working directory
#output_dir = os.path.join(os.getcwd(), "assembly_output")


# Define the MegaHit command
#megahit_cmd = f"megahit -1 {','.join(R1_files)} -2 {','.join(R2_files)} -o {output_dir} --k-min 21 --k-max 101 --min-count 2"

# Run MegaHit using subprocess
#subprocess.run(megahit_cmd, shell=True)

######**************** generate gfa graph from final config fa file from megahit

#fa2gfa -i final.contigs.fa -o final_contigs.gfa



#import csv
#.........................................................................................................
#def extract_contig_info(gfa_file, csv_file):
   # with open(gfa_file, 'r') as gfa, open(csv_file, 'w', newline='') as csv_out:
       # csv_writer = csv.writer(csv_out)
       # csv_writer.writerow(['name', 'length'])

       # for line in gfa:
           # if line.startswith('S'):
              #  parts = line.strip().split('\t')
               # name = parts[1]
               # sequence = parts[2]
              #  length = len(sequence)
               # csv_writer.writerow([name, length])

# Example usage
# gfa_file = 'final_contigs.gfa'
#csv_file = 'annotations.csv'
#extract_contig_info(gfa_file, csv_file)
#print(f"CSV file '{csv_file}' created successfully.")

#..........................................................................................................
#csv2tag -i final_contigs.gfa -c annotations.csv -o tagged.gfa


# ************************** script to generate gfa with nodes connected to edges
#def add_links_to_gfa(input_gfa, output_gfa):
   # segments = []
   # with open(input_gfa, 'r') as gfa, open(output_gfa, 'w') as out_gfa:
        #for line in gfa:
           # parts = line.strip().split('\t')
            #if parts[0] == 'S':
               # segments.append(parts[1])
            #out_gfa.write(line)

        ########## Add example links
       # for i in range(len(segments) - 1):
          #  out_gfa.write(f"L\t{segments[i]}\t+\t{segments[i+1]}\t+\t0M\n")

# Example usage
#input_gfa = 'final_contigs.gfa'
#output_gfa = 'final_contigs_with_links.gfa'
#add_links_to_gfa(input_gfa, output_gfa)
#print(f"Updated GFA file with links saved as '{output_gfa}'.")

















