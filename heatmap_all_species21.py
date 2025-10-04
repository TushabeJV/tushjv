import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the data
filename = 'kraken2_report_with_rank_names.csv'
df = pd.read_csv(filename)

# Filter for species level (taxon_rank 'S') and get species with ≥50 unique reads
species_df = df[df['taxon_rank'] == 'S'].copy()
filtered_species = species_df[species_df['unique_read_count'] >= 50]

print(f"Total species with ≥50 unique reads: {len(filtered_species)}")

if len(filtered_species) == 0:
    print("No species found with ≥50 unique reads!")
else:
    # Sort by unique read count for better visualization
    filtered_species = filtered_species.sort_values('unique_read_count', ascending=False)

    # OPTION 1: Create a paginated version for large datasets
    def create_paginated_heatmap(data, species_per_page=50):
        num_pages = (len(data) + species_per_page - 1) // species_per_page

        for page in range(num_pages):
            start_idx = page * species_per_page
            end_idx = min((page + 1) * species_per_page, len(data))
            page_data = data.iloc[start_idx:end_idx]

            fig, ax = plt.subplots(figsize=(16, min(12, species_per_page * 0.4)))

            # Create horizontal bar plot as heatmap
            norm = plt.Normalize(page_data['unique_read_count'].min(),
                                 page_data['unique_read_count'].max())
            colors = plt.cm.viridis(norm(page_data['unique_read_count']))

            bars = ax.barh(range(len(page_data)),
                           page_data['unique_read_count'],
                           color=colors,
                           edgecolor='black',
                           linewidth=0.5)

            # Customize the plot
            ax.set_title(f'Species by Unique Read Count (Page {page+1}/{num_pages})\nSpecies {start_idx+1}-{end_idx} of {len(data)}',
                         fontsize=14, fontweight='bold', pad=20)
            ax.set_xlabel('Unique Read Count', fontsize=12)
            ax.set_ylabel('Species', fontsize=12)

            # Set y-axis ticks to species names
            ax.set_yticks(range(len(page_data)))
            ax.set_yticklabels(page_data['name'], fontsize=8)

            # Add value labels on the right
            max_reads = page_data['unique_read_count'].max()
            for i, (bar, count) in enumerate(zip(bars, page_data['unique_read_count'])):
                ax.text(bar.get_width() + max_reads * 0.01,
                        bar.get_y() + bar.get_height()/2,
                        f"{count:,}",
                        va='center', ha='left', fontsize=7)

            # Add rank code at the bottom
            ax.text(0.5, -0.15, 'Rank Code: S', transform=ax.transAxes,
                    fontsize=10, ha='center', va='top')

            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.8, aspect=20, pad=0.02)
            cbar.set_label('Unique Read Count', rotation=270, labelpad=15)

            # Use log scale if data range is large
            if page_data['unique_read_count'].max() / page_data['unique_read_count'].min() > 100:
                ax.set_xscale('log')
                ax.set_xlabel('Unique Read Count (Log Scale)', fontsize=12)

            plt.tight_layout()

            output_filename = f'species_heatmap_page_{page+1}.png'
            plt.savefig(output_filename, dpi=200, bbox_inches='tight')
            print(f"Page {page+1} saved as: {output_filename}")
            plt.close(fig)  # Closing the figure to free memory

    # Create paginated heatmaps
    create_paginated_heatmap(filtered_species, species_per_page=50)

    # OPTION 2: Create a condensed version with top N species
    print("\nCreating condensed version with top species...")

    # Show top 50 species in one plot
    top_species = filtered_species.head(50)

    fig, ax = plt.subplots(figsize=(14, 12))

    norm = plt.Normalize(top_species['unique_read_count'].min(),
                         top_species['unique_read_count'].max())
    colors = plt.cm.plasma(norm(top_species['unique_read_count']))

    bars = ax.barh(range(len(top_species)),
                   top_species['unique_read_count'],
                   color=colors,
                   edgecolor='black',
                   linewidth=0.5)

    ax.set_title('Top 50 Species by Unique Read Count',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Unique Read Count', fontsize=12)
    ax.set_ylabel('Species', fontsize=12)

    ax.set_yticks(range(len(top_species)))
    ax.set_yticklabels(top_species['name'], fontsize=9)

    # Add value labels
    max_reads = top_species['unique_read_count'].max()
    for i, (bar, count) in enumerate(zip(bars, top_species['unique_read_count'])):
        ax.text(bar.get_width() + max_reads * 0.01,
                bar.get_y() + bar.get_height()/2,
                f"{count:,}",
                va='center', ha='left', fontsize=8)

    ax.text(0.5, -0.05, 'Rank Code: S', transform=ax.transAxes,
            fontsize=12, ha='center', va='top')

    sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Unique Read Count', rotation=270, labelpad=15)

    if top_species['unique_read_count'].max() / top_species['unique_read_count'].min() > 100:
        ax.set_xscale('log')
        ax.set_xlabel('Unique Read Count (Log Scale)', fontsize=12)

    plt.tight_layout()
    plt.savefig('species_heatmap_top50.png', dpi=300, bbox_inches='tight')
    print("Top 50 heatmap saved as: species_heatmap_top50.png")
    plt.show()

    # OPTION 3: Creating an interactive HTML version of all data
    try:
        import plotly.express as px
        import plotly.graph_objects as go

        # Creating interactive horizontal bar chart
        fig = px.bar(filtered_species,
                     x='unique_read_count',
                     y='name',
                     orientation='h',
                     title='Classified Species',
                     labels={'unique_read_count': 'Unique Read Count', 'name': 'Species'},
                     color='unique_read_count',
                     color_continuous_scale='viridis')

        fig.update_layout(
            height=max(600, len(filtered_species) * 20),
            yaxis={'categoryorder': 'total ascending'},
            showlegend=False
        )

        fig.write_html('species_interactive_heatmap.html')
        print("Interactive HTML version saved as: species_interactive_heatmap.html")

    except ImportError:
        print("Plotly not available for interactive version. Install with: pip install plotly")

    # Save detailed results
    results_df = filtered_species[['name', 'unique_read_count', 'read_count', 'percentage_abundance']].copy()
    results_df = results_df.sort_values('unique_read_count', ascending=False)
    results_df.columns = ['Species', 'Unique_Reads', 'Total_Reads', 'Percentage_Abundance']
    results_df.to_csv('species_50plus_unique_reads_detailed.csv', index=False)
    print(f"\nDetailed results saved to: species_50plus_unique_reads_detailed.csv")

    # Printing summary
    print(f"\nDataset summary:")
    print(f"  - Total species with ≥50 unique reads: {len(filtered_species)}")
    print(f"  - Total unique reads in filtered species: {filtered_species['unique_read_count'].sum():,}")
    print(f"  - Range: {filtered_species['unique_read_count'].min():,} to {filtered_species['unique_read_count'].max():,}")