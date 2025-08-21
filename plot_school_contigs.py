import matplotlib.pyplot as plt
import re
import numpy as np

# Initialize lists to store locations and total contigs
locations = []
total_contigs = []

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
            # Remove commas and convert to integer
            contigs = int(contigs_match.group(1).replace(',', ''))
            total_contigs.append(contigs)

# Generate unique colors using a colormap
num_schools = len(locations)  # Expecting 21 schools
colors = plt.cm.tab20(np.linspace(0, 1, min(num_schools, 20)))  # tab20 colormap for up to 20 colors
if num_schools > 20:
    # Extend with additional distinct colors if more than 20 schools
    extra_colors = plt.cm.Set2(np.linspace(0, 1, num_schools - 20))
    colors = np.vstack([colors, extra_colors])

# Calculate the average of total contigs
average_contigs = np.mean(total_contigs)

# Create the bar plot
plt.figure(figsize=(12, 6))  # Increased width to fit 21 bars
bars = plt.bar(locations, total_contigs, color=colors)
plt.title("Total Contigs by School")
plt.xlabel("School")
plt.ylabel("Total Contigs")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for readability

# Add horizontal line for the average
plt.axhline(y=average_contigs, color='black', linestyle='--', linewidth=2, label=f'Average: {average_contigs:,.0f}')

# Add text label on the line
plt.text(x=len(locations)-0.5, y=average_contigs + (max(total_contigs) - min(total_contigs)) * 0.02,
         s='Average contigs', ha='right', va='bottom', fontsize=10, color='black')

# Add a legend
plt.legend(bars, locations + [f'Average: {average_contigs:,.0f}'], title="Schools", bbox_to_anchor=(1.05, 1), loc='upper left')

# Display the plot
plt.tight_layout()
plt.show()