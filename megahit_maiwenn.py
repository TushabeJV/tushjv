import subprocess
import os
import shutil

nanopore_fastq = "all_reads.fastq"
output_dir = os.path.join(os.getcwd(), "megahit_nanopore_output")

# remove the existing output directory, if any
if os.path.exists(output_dir):
    print(f"Removing existing output directory: {output_dir}")
    shutil.rmtree(output_dir)  # Important: let MEGAHIT recreate it


megahit_cmd = [
    "megahit",
    "-r", nanopore_fastq,
    "-o", output_dir,
    "--min-count", "2",
    "--k-min", "21",
    "--k-max", "101",
    "--k-step", "10",
    "--presets", "meta-sensitive"
]

print(f"Running MEGAHIT on {nanopore_fastq}...")
result = subprocess.run(megahit_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

if result.returncode != 0:
    print("MEGAHIT failed:")
    print(result.stderr)
else:
    print("Assembly complete.")
    print(f"Contigs saved in: {os.path.join(output_dir, 'final.contigs.fa')}")
