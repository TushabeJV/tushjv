def ensure_links_in_gfa(input_gfa, output_gfa):
    segments = []
    with open(input_gfa, 'r') as gfa, open(output_gfa, 'w') as out_gfa:
        for line in gfa:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                segments.append(parts[1])
            out_gfa.write(line)

        # Add example links between segments
        for i in range(len(segments) - 1):
            out_gfa.write(f"L\t{segments[i]}\t+\t{segments[i+1]}\t+\t0M\n")

# Example usage
input_gfa = 'final_contigs.gfa'
output_gfa = 'final_contigs_with_links2.gfa'
ensure_links_in_gfa(input_gfa, output_gfa)
print(f"Updated GFA file with links saved as '{output_gfa}'.")
