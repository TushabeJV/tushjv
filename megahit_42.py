import os
import glob
import subprocess

# Configuration
input_dir = "Trimmed_42"
individual_out_dir = "megahit_individual"
combined_out_dir = "megahit_combined"
megahit_params = "--k-min 21 --k-max 101 --min-count 2"

# Create output directories
os.makedirs(individual_out_dir, exist_ok=True)
os.makedirs(combined_out_dir, exist_ok=True)

# Get all R1 files and extract unique samples
r1_files = sorted(glob.glob(os.path.join(input_dir, "*_R1_*.fq.gz")))
samples = set()
for r1_file in r1_files:
    filename = os.path.basename(r1_file)
    sample_id = "_".join(filename.split("_")[:4])  # e.g., "N01_S_2_S18"
    samples.add(sample_id)

print(f"Found {len(samples)} unique samples for individual assembly")
print(f"Found {len(r1_files)} R1 files for co-assembly")

# 1. Run individual assemblies
for sample in sorted(samples):
    print(f"\nStarting assembly for sample: {sample}")

    # Find sample-specific files
    sample_r1 = glob.glob(os.path.join(input_dir, f"{sample}_R1_*.fq.gz"))
    sample_r2 = glob.glob(os.path.join(input_dir, f"{sample}_R2_*.fq.gz"))

    if not sample_r1 or not sample_r2:
        print(f"  ! Missing files for sample {sample}, skipping")
        continue

    # Build MEGAHIT command
    cmd = (
        f"megahit -1 {','.join(sample_r1)} -2 {','.join(sample_r2)} "
        f"-o {os.path.join(individual_out_dir, sample)} "
        f"{megahit_params}"
    )

    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"  ✓ Completed {sample} assembly")
    else:
        print(f"  ! Failed {sample} assembly. Error:\n{result.stderr}")

# 2. Run co-assembly of all samples
print("\nStarting co-assembly of all samples")
all_r1 = glob.glob(os.path.join(input_dir, "*_R1_*.fq.gz"))
all_r2 = glob.glob(os.path.join(input_dir, "*_R2_*.fq.gz"))

cmd = (
    f"megahit -1 {','.join(all_r1)} -2 {','.join(all_r2)} "
    f"-o {combined_out_dir} {megahit_params}"
)

print(f"Running: {cmd}")
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print("✓ Completed co-assembly")
    # Check contig output
    contig_file = os.path.join(combined_out_dir, "final.contigs.fa")
    if os.path.exists(contig_file):
        print(f"  Contigs generated: {contig_file}")
    else:
        print("  ! Contig file missing - check MEGAHIT logs")
else:
    print(f"! Co-assembly failed. Error:\n{result.stderr}")

print("\nAll assemblies completed!")