import os
import glob
import subprocess
from collections import defaultdict

def find_fastq_pairs(fastq_dir):
    """Automatically detect and pair R1 and R2 FASTQ files based on sample ID."""
    fastq_files = glob.glob(os.path.join(fastq_dir, "*_R[12]_001.fastq.gz"))

    # Group files by sample ID (everything before _R1 or _R2)
    sample_dict = defaultdict(dict)

    for file in fastq_files:
        basename = os.path.basename(file)
        # Extract sample ID (everything before _R1_ or _R2_)
        if "_R1_" in basename:
            sample_id = basename.split("_R1_")[0]
            sample_dict[sample_id]["R1"] = file
        elif "_R2_" in basename:
            sample_id = basename.split("_R2_")[0]
            sample_dict[sample_id]["R2"] = file

    # Verify we have both R1 and R2 for each sample
    forward_reads = []
    reverse_reads = []
    unmatched_samples = []

    for sample_id, files in sample_dict.items():
        if "R1" in files and "R2" in files:
            forward_reads.append(files["R1"])
            reverse_reads.append(files["R2"])
        else:
            unmatched_samples.append(sample_id)

    if unmatched_samples:
        print(f"Warning: {len(unmatched_samples)} samples missing R1/R2 pairs:")
        for sample in unmatched_samples[:5]:  # Print first 5 to avoid flooding
            print(f"  {sample}")
        if len(unmatched_samples) > 5:
            print(f"  ...and {len(unmatched_samples)-5} more")

    return forward_reads, reverse_reads

def run_trimmomatic(forward_reads, reverse_reads, output_dir, adapter_file, quality="phred33"):
    """Run Trimmomatic on detected FASTQ pairs with enhanced error handling."""
    os.makedirs(output_dir, exist_ok=True)

    # Verify Trimmomatic is available
    try:
        subprocess.run(["trimmomatic", "-version"], check=True, capture_output=True)
    except FileNotFoundError:
        raise RuntimeError("Trimmomatic not found. Please ensure it's installed and in your PATH.")

    for fwd, rev in zip(forward_reads, reverse_reads):
        sample_id = os.path.basename(fwd).split("_R1_")[0]

        # Output file names
        out_files = [
            os.path.join(output_dir, f"{sample_id}_R1_trimmed_paired.fq.gz"),
            os.path.join(output_dir, f"{sample_id}_R1_trimmed_unpaired.fq.gz"),
            os.path.join(output_dir, f"{sample_id}_R2_trimmed_paired.fq.gz"),
            os.path.join(output_dir, f"{sample_id}_R2_trimmed_unpaired.fq.gz")
        ]

        # Skip if output files already exist
        if all(os.path.exists(f) for f in out_files):
            print(f"✓ Output files exist, skipping {sample_id}")
            continue

        # Trimmomatic command with enhanced parameters
        cmd = [
            "trimmomatic", "PE",
            "-threads", str(os.cpu_count()),
            f"-{quality}",
            fwd, rev,
            *out_files,
            f"ILLUMINACLIP:{adapter_file}:2:30:10:8:true",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36",
            "-validatePairs"  # Ensure paired reads match
        ]

        print(f"\nProcessing {sample_id}...")
        print(f"  Input: {os.path.basename(fwd)} | {os.path.basename(rev)}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✓ Successfully processed {sample_id}")

            # Log trimming results
            for line in result.stdout.split('\n'):
                if "Input Read Pairs" in line:
                    print(f"  {line.strip()}")
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to process {sample_id}")
            print(f"  Error: {e.stderr.strip()}")
            # Clean up partial outputs if any
            for f in out_files:
                if os.path.exists(f):
                    os.remove(f)

if __name__ == "__main__":
    # Configuration - modify these as needed
    CONFIG = {
        "fastq_dir": "Fastq_42",
        "output_dir": "Trimmed_Reads",
        "adapter_file": "NexteraPE-PE.fa",
        "quality": "phred33",
        "min_file_size": 20          # Minimum file size in bytes (to skip empty files)
    }

    print("🔍 Detecting FASTQ file pairs...")
    try:
        forward_reads, reverse_reads = find_fastq_pairs(CONFIG["fastq_dir"])
        print(f"✔ Found {len(forward_reads)} valid sample pairs.")

        # Filter out suspiciously small files
        filtered_pairs = []
        for fwd, rev in zip(forward_reads, reverse_reads):
            fwd_size = os.path.getsize(fwd)
            rev_size = os.path.getsize(rev)
            if fwd_size < CONFIG["min_file_size"] or rev_size < CONFIG["min_file_size"]:
                sample_id = os.path.basename(fwd).split("_R1_")[0]
                print(f"⚠ Skipping {sample_id} - file too small (FWD: {fwd_size}b, REV: {rev_size}b)")
                continue
            filtered_pairs.append((fwd, rev))

        if filtered_pairs:
            forward_reads, reverse_reads = zip(*filtered_pairs)
            print(f"Processing {len(forward_reads)} samples after size filtering...")
            run_trimmomatic(
                forward_reads=forward_reads,
                reverse_reads=reverse_reads,
                output_dir=CONFIG["output_dir"],
                adapter_file=CONFIG["adapter_file"],
                quality=CONFIG["quality"]
            )
        else:
            print("No valid sample pairs found after filtering.")

    except Exception as e:
        print(f"❌ Error: {str(e)}")
    finally:
        print("\nProcessing complete.")