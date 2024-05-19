# viewing fastq stats
import subprocess
import os
#defining input file in output directory

#input_files = ["Fastq/BO4_S1_L001_R1_001.fastq", "Fastq/BO4_S1_L001_R2_001.fastq", "Fastq/B05_S0_L001_R1_001.fastq", "Fastq/BO5_S0_L001_R2_001.fastq",
#"Fastq/G03_S3_L001_R1_001.fastq", "Fastq/G03_S3_L001_R2_001.fastq", "Fastq/G11_S4_L001_R1_001.fastq", "Fastq/G11_S4_L001_R2_001.fastq" ]

input_files = ["BO4_trimmed_forward_paired.fq.gz", "BO4_trimmed_reverse_paired.fq.gz" ]
output_dir = "Fastq"

# run fastqc using suprocess
for input_file in input_files:
 subprocess.run(f"fastqc {input_file} -o {output_dir}", shell=True)


























