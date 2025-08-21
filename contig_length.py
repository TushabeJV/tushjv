import matplotlib.pyplot as plt
import re
import numpy as np

# Initialize lists to store locations, total contigs, and total lengths
locations = []
total_contigs = []
total_lengths = []

# Read the summary_report.txt file
with open('summary_report.txt', 'r') as file:
    lines = file.readlines()
    for line in lines:
        # Match lines with "Location" to get the location name
        location_match = re.match(r'Location (\w+):', line.strip())
        if location_match:
            current_location = location_match.group(1)
            locations.append(current_location)
        # Match lines with "Total contigs" to get the contig count
        contigs_match = re.match(r'\s*Total contigs: ([\d,]+)', line.strip())
        if contigs_match:
            contigs = int(contigs_match.group(1).replace(',', ''))
            total_contigs.append(contigs)
        # Match lines with "Total length" to get the length in bp
        length_match = re.match(r'\s*Total length: ([\d,]+) bp', line.strip())
        if length_match:
            length = int(length_match.group(1).replace(',', ''))
            total_lengths.append(length)

# Single colors for bars
contigs_color = '#2a623d'  # Dark green for all Total Contigs bars
length_color = '#1f4e79'   # Dark blue for all Total Length bars

# Calculate averages
average_contigs = np.mean(total_contigs)
average_length = np.mean(total_lengths)

# Create the figure and dual axes
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot Total Contigs on the primary y-axis (left)
bars_contigs = ax1.bar(np.arange(len(locations)) - 0.2, total_contigs, width=0.4, color=contigs_color, label='Total Contigs')
ax1.set_xlabel('School')
ax1.set_ylabel('Total Contigs', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xticks(np.arange(len(locations)))
ax1.set_xticklabels(locations, rotation=45, ha='right')
ax1.grid(axis='y', linestyle='--', alpha=0.7)

# Add horizontal line and label for average contigs (left side)
ax1.axhline(y=average_contigs, color='black', linestyle='--', linewidth=2)
ax1.text(x=0.5, y=average_contigs + (max(total_contigs) - min(total_contigs)) * 0.03,
         s='Average contigs', ha='left', va='bottom', fontsize=10, color='black')

# Create secondary y-axis for Total Length
ax2 = ax1.twinx()
bars_length = ax2.bar(np.arange(len(locations)) + 0.2, total_lengths, width=0.4, color=length_color, alpha=0.5, label='Total Length')
ax2.set_ylabel('Total Length (bp)', color='darkblue')
ax2.tick_params(axis='y', labelcolor='darkblue')

# Add horizontal line and label for average length (right side)
ax2.axhline(y=average_length, color='darkred', linestyle='--', linewidth=2)
ax2.text(x=len(locations)-0.5, y=average_length + (max(total_lengths) - min(total_lengths)) * 0.03,
         s='Average length', ha='right', va='bottom', fontsize=10, color='darkred')

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2 + [plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2),
                              plt.Line2D([0], [0], color='darkred', linestyle='--', linewidth=2)],
           labels1 + labels2 + [f'Average Contigs: {average_contigs:,.0f}', f'Average Length: {average_length:,.0f}'],
           title="Metrics", bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust layout and display
plt.title("Total Contigs and Total Length by School")
plt.tight_layout()
plt.show()