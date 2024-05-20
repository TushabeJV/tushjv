import subprocess
import os

# Define the input files
R1_files = ["Trimmed_paired/BO4_trimmed_r1_paired.fq.gz", "Trimmed_paired/B05_trimmed_r1_paired.fq.gz", "Trimmed_paired/G03_trimmed_r1_paired.fq.gz", "Trimmed_paired/G11_trimmed_r1_paired.fq.gz"]
R2_files = ["Trimmed_paired/BO4_trimmed_r2_paired.fq.gz", "Trimmed_paired/B05_trimmed_r2_paired.fq.gz", "Trimmed_paired/G03_trimmed_r2_paired.fq.gz", "Trimmed_paired/G11_trimmed_r2_paired.fq.gz"]

# Set the output directory to a new folder called "assembly_output" in the current working directory
output_dir = os.path.join(os.getcwd(), "assembly_output")


# Define the MegaHit command
megahit_cmd = f"megahit -1 {','.join(R1_files)} -2 {','.join(R2_files)} -o {output_dir} --k-min 21 --k-max 101 --min-count 2"

# Run MegaHit using subprocess
subprocess.run(megahit_cmd, shell=True)
