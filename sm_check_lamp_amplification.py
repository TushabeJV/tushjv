from Bio import SeqIO

# Load your target sequence (FASTA format)
target_file = "final_contigs.fa"  # Your target sequence file
primers = {
    'F3': "TCGTCTATAGTACGGTAGG",
    'B3': "TTATACTTTAACCCCCACC",
    'FIP': "TGCCAAGTAGAGACTACAAACATCTTGGGTAAGGTAGAAAATGTT",
    'BIP': "AGAAGTGTTTAACTTGATGAAGGGGAAACAAAACCGAAACCA",
    'LF': "GTCCTCTTGTTTTTGAAT",
    'LB': "CTGCACGAAATACAGAATC"
}

# Reading the target sequences
for record in SeqIO.parse(target_file, "fasta"):
    target_seq = str(record.seq)

    # Check if primers are found in the target sequence
    primer_matches = {}

    # Check F3
    if primers['F3'] in target_seq:
        primer_matches['F3'] = "Forward primer (F3) binds"

    # Check B3
    if primers['B3'] in target_seq:
        primer_matches['B3'] = "Backward primer (B3) binds"

    # Check FIP
    if primers['FIP'] in target_seq:
        primer_matches['FIP'] = "Forward internal primer (FIP) binds"

    # Check BIP
    if primers['BIP'] in target_seq:
        primer_matches['BIP'] = "Backward internal primer (BIP) binds"

    # Check LF
    if primers['LF'] in target_seq:
        primer_matches['LF'] = "Loop forward primer (LF) binds"

    # Check LB
    if primers['LB'] in target_seq:
        primer_matches['LB'] = "Loop backward primer (LB) binds"

    # Print results only if any primer binds
    if primer_matches:
        print(f"Target sequence: {record.id}")
        for primer, result in primer_matches.items():
            print(f"{primer}: {result}")
