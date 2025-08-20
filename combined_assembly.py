import os
import glob
import subprocess
import shutil
import time

# Configuration
input_dir = "Trimmed_Reads/Trimmed_42"
combined_out_dir = "megahit_combined"
megahit_params = "--k-min 21 --k-max 101 --min-count 2"

# Function to safely remove directory with retries
def safe_remove_directory(dir_path, max_retries=5, delay=2):
    for attempt in range(max_retries):
        try:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                print(f"Successfully removed directory: {dir_path}")
                return True
            return True
        except Exception as e:
            print(f"Attempt {attempt + 1} failed to remove {dir_path}: {str(e)}")
            if attempt < max_retries - 1:
                print(f"Waiting {delay} seconds before retry...")
                time.sleep(delay)
    return False

# Try to remove existing directory
if not safe_remove_directory(combined_out_dir):
    print(f"Could not remove {combined_out_dir}, using alternative approach")
    # Create a unique directory name with timestamp
    import datetime
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    combined_out_dir = f"megahit_combined_{timestamp}"

# Get all R1 and R2 files for co-assembly
all_r1 = sorted(glob.glob(os.path.join(input_dir, "*_R1_*.fq.gz")))
all_r2 = sorted(glob.glob(os.path.join(input_dir, "*_R2_*.fq.gz")))

print(f"Found {len(all_r1)} R1 files and {len(all_r2)} R2 files for co-assembly")

# Verify we have matching numbers of R1 and R2 files
if len(all_r1) != len(all_r2):
    print("Warning: Number of R1 and R2 files don't match!")
    print("This might indicate missing or misnamed files")

# Run co-assembly of all samples
print("\nStarting co-assembly of all samples")

cmd = (
    f"megahit -1 {','.join(all_r1)} -2 {','.join(all_r2)} "
    f"-o {combined_out_dir} {megahit_params}"
)

print(f"Running: {cmd}")
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print("<E2><9C><93> Completed co-assembly")
    # Check contig output
    contig_file = os.path.join(combined_out_dir, "final.contigs.fa")
    if os.path.exists(contig_file):
        contig_size = os.path.getsize(contig_file)
        print(f"  Contigs generated: {contig_file} ({contig_size} bytes)")

        # Count the number of contigs
        try:
            contig_count = subprocess.check_output(f"grep -c '>' {contig_file}", shell=True).decode().strip()
            print(f"  Number of contigs: {contig_count}")
        except:
            print("  Could not count contigs")
    else:
        print("  ! Contig file missing - check MEGAHIT logs")
        if result.stdout:
            print(f"  MEGAHIT output: {result.stdout[-500:]}")
else:
    print(f"! Co-assembly failed. Error:\n{result.stderr}")
    if result.stdout:
        print(f"MEGAHIT output: {result.stdout[-500:]}")

print("\nCo-assembly completed!")