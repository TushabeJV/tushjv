import os
import subprocess

def run_trimmomatic(forward_reads, reverse_reads, output_dir, adapter_file):
    # defining output_dir, works better than just defining directory
    if not output_dir:
        output_dir = os.getcwd()

    # Iterate over each pair of forward and reverse reads
    for i, (forward_read, reverse_read) in enumerate(zip(forward_reads, reverse_reads)):

        # Extract sample identifier from the forward read filename
        sample_id = os.path.basename(forward_read).split('_')[0]


        # Create output file names
        forward_paired = os.path.join(output_dir, f"{sample_id}_trimmed_forward_paired.fq.gz")
        reverse_paired = os.path.join(output_dir, f"{sample_id}_trimmed_reverse_paired.fq.gz")
        forward_unpaired = os.path.join(output_dir, f"{sample_id}_trimmed_forward_unpaired.fq.gz")
        reverse_unpaired = os.path.join(output_dir, f"{sample_id}_trimmed_reverse_unpaired.fq.gz")

        # Construct the Trimmomatic command
        trimmomatic_command = ["trimmomatic", "PE", "-phred33", forward_read, reverse_read, forward_paired, forward_unpaired, reverse_paired, reverse_unpaired, f"ILLUMINACLIP:{adapter_file}:2:30:10",
            "SLIDINGWINDOW:4:15", "MINLEN:36"]

        # Print the command for debugging
        print(f"Running command for sample {sample_id}:", " ".join(trimmomatic_command))

        # Execute the command using subprocess
        try:
            result = subprocess.run(trimmomatic_command, check=True, capture_output=True, text=True)
            print(f"Trimmomatic output for sample {sample_id}:", result.stdout)
            print(f"Trimmomatic errors for sample {sample_id}:", result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error running Trimmomatic for sample {sample_id}:", e.stderr)

# Example usage
forward_reads = ["Fastq/BO4_S1_R1.fastq.gz", "Fastq/B05_S0_R1.fastq.gz", "Fastq/G03_S3_R1.fastq.gz", "Fastq/G11_S4_R1.fastq.gz"]
reverse_reads = ["Fastq/BO4_S1_R2.fastq.gz", "Fastq/B05_S0_R2.fastq.gz", "Fastq/G03_S3_R2.fastq.gz", "Fastq/G11_S4_R2.fastq.gz"]

run_trimmomatic(forward_reads=forward_reads, reverse_reads=reverse_reads, output_dir=os.getcwd(), adapter_file="TruSeq3_PE.fa")

