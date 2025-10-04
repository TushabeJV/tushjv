import pandas as pd
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
    # Sort by unique read count for better visualization (ascending for horizontal bars)
    filtered_species = filtered_species.sort_values('unique_read_count', ascending=True)

    try:
        import plotly.express as px
        import plotly.graph_objects as go

        # Create interactive horizontal bar chart with reduced bar width
        fig = px.bar(filtered_species,
                     x='unique_read_count',
                     y='name',
                     orientation='h',
                     title='Classified Species',
                     labels={'unique_read_count': 'Unique Read Count', 'name': 'Species'},
                     color='unique_read_count',
                     color_continuous_scale='viridis',
                     hover_data={
                         'name': True,
                         'unique_read_count': True,
                         'read_count': True,
                         'percentage_abundance': ':.4f'
                     })

        # Update layout for better appearance with thinner bars
        fig.update_layout(
            height=max(800, len(filtered_species) * 15),  # Reduced from 25 to 15 for thinner bars
            width=1200,
            yaxis={
                'categoryorder': 'total ascending',
                'title': 'Species',
                'tickfont': {'size': 9},
                'automargin': True
            },
            xaxis={
                'title': 'Unique Read Count',
                'tickformat': ','
            },
            showlegend=False,
            title_x=0.5,  # Center the title
            title_font_size=24,
            title_font_color='darkblue',
            font=dict(family="Arial, sans-serif", size=11),
            plot_bgcolor='rgba(248,248,248,1)',
            paper_bgcolor='rgba(248,248,248,1)',
            bargap=0.1,  # Reduced gap between bars
        )

        # Update traces to make bars thinner
        fig.update_traces(
            marker=dict(
                line=dict(width=0.5, color='DarkSlateGrey')
            ),
            width=0.6,  # Reduced bar width (default is 0.8)
            hovertemplate='<b>%{y}</b><br>' +
                          'Unique Reads: %{x:,}<br>' +
                          'Total Reads: %{customdata[0]:,}<br>' +
                          'Percentage: %{customdata[1]:.4f}%<br>' +
                          '<extra></extra>'
        )

        # Add a subtitle
        fig.add_annotation(
            text=f"{len(filtered_species)} Species with ≥50 Unique Reads",
            xref="paper", yref="paper",
            x=0.5, y=1.02,
            xanchor='center',
            showarrow=False,
            font=dict(size=14, color="gray")
        )

        # Add rank code annotation at the bottom
        fig.add_annotation(
            text="Rank Code: S",
            xref="paper", yref="paper",
            x=0.5, y=-0.08,
            xanchor='center',
            showarrow=False,
            font=dict(size=12, color="black")
        )

        # Save the interactive plot
        output_filename = 'classified_species_thin_bars.html'
        fig.write_html(output_filename)
        print(f"Interactive plot with thin bars saved as: {output_filename}")

        # Create an alternative version with even thinner bars
        fig_thinner = px.bar(filtered_species,
                             x='unique_read_count',
                             y='name',
                             orientation='h',
                             title='',
                             labels={'unique_read_count': 'Unique Read Count', 'name': 'Species'},
                             color='unique_read_count',
                             color_continuous_scale='plasma',
                             hover_data={
                                 'name': True,
                                 'unique_read_count': True,
                                 'read_count': True,
                                 'percentage_abundance': ':.4f'
                             })

        # Even thinner bars
        fig_thinner.update_layout(
            height=max(800, len(filtered_species) * 12),  # Even more compact
            width=1200,
            yaxis={
                'categoryorder': 'total ascending',
                'title': 'Species',
                'tickfont': {'size': 8},
                'automargin': True
            },
            xaxis={
                'title': 'Unique Read Count',
                'tickformat': ','
            },
            showlegend=False,
            title_x=0.5,
            title_font_size=22,
            title_font_color='darkred',
            font=dict(family="Arial, sans-serif", size=10),
            plot_bgcolor='rgba(248,248,248,1)',
            paper_bgcolor='rgba(248,248,248,1)',
            bargap=0.15,  # Slightly more gap for very thin bars
        )

        fig_thinner.update_traces(
            marker=dict(
                line=dict(width=0.3, color='DarkSlateGrey')
            ),
            width=0.4,  # Even thinner bars
            hovertemplate='<b>%{y}</b><br>' +
                          'Unique Reads: %{x:,}<br>' +
                          'Total Reads: %{customdata[0]:,}<br>' +
                          'Percentage: %{customdata[1]:.4f}%<br>' +
                          '<extra></extra>'
        )

        fig_thinner.add_annotation(
            text=f"{len(filtered_species)} Species with ≥50 Unique Reads ",
            xref="paper", yref="paper",
            x=0.5, y=1.02,
            xanchor='center',
            showarrow=False,
            font=dict(size=12, color="gray")
        )

        fig_thinner.add_annotation(
            text="Rank Code: S",
            xref="paper", yref="paper",
            x=0.5, y=-0.08,
            xanchor='center',
            showarrow=False,
            font=dict(size=11, color="black")
        )

        thinner_output_filename = 'classified_species_very_thin_bars.html'
        fig_thinner.write_html(thinner_output_filename)
        print(f"Interactive plot with very thin bars saved as: {thinner_output_filename}")

        # Create a log scale version with thin bars
        fig_log = px.bar(filtered_species,
                         x='unique_read_count',
                         y='name',
                         orientation='h',
                         title='Classified Species (Log Scale)',
                         labels={'unique_read_count': 'Unique Read Count (Log Scale)', 'name': 'Species'},
                         color='unique_read_count',
                         color_continuous_scale='viridis',
                         hover_data={
                             'name': True,
                             'unique_read_count': True,
                             'read_count': True,
                             'percentage_abundance': ':.4f'
                         })

        fig_log.update_layout(
            height=max(800, len(filtered_species) * 15),
            width=1200,
            yaxis={
                'categoryorder': 'total ascending',
                'title': 'Species',
                'tickfont': {'size': 9},
                'automargin': True
            },
            xaxis={
                'title': 'Unique Read Count (Log Scale)',
                'type': 'log'
            },
            showlegend=False,
            title_x=0.5,
            title_font_size=22,
            font=dict(family="Arial, sans-serif", size=11),
            plot_bgcolor='rgba(248,248,248,1)',
            paper_bgcolor='rgba(248,248,248,1)',
            bargap=0.1,
        )

        fig_log.update_traces(
            marker=dict(
                line=dict(width=0.5, color='DarkSlateGrey')
            ),
            width=0.6,
            hovertemplate='<b>%{y}</b><br>' +
                          'Unique Reads: %{x:,}<br>' +
                          'Total Reads: %{customdata[0]:,}<br>' +
                          'Percentage: %{customdata[1]:.4f}%<br>' +
                          '<extra></extra>'
        )

        fig_log.add_annotation(
            text=f"{len(filtered_species)} Species with ≥50 Unique Reads",
            xref="paper", yref="paper",
            x=0.5, y=1.02,
            xanchor='center',
            showarrow=False,
            font=dict(size=14, color="gray")
        )

        fig_log.add_annotation(
            text="Rank Code: S",
            xref="paper", yref="paper",
            x=0.5, y=-0.08,
            xanchor='center',
            showarrow=False,
            font=dict(size=12, color="black")
        )

        log_output_filename = 'classified_species_thin_bars_log.html'
        fig_log.write_html(log_output_filename)
        print(f"Interactive plot with thin bars (log scale) saved as: {log_output_filename}")

        # Create a summary statistics table
        print(f"\n{'='*60}")
        print("DATASET SUMMARY")
        print(f"{'='*60}")
        print(f"Total species with ≥50 unique reads: {len(filtered_species):,}")
        print(f"Total unique reads in dataset: {filtered_species['unique_read_count'].sum():,}")
        print(f"Average unique reads per species: {filtered_species['unique_read_count'].mean():,.1f}")
        print(f"Median unique reads per species: {filtered_species['unique_read_count'].median():,.1f}")
        print(f"Range: {filtered_species['unique_read_count'].min():,} - {filtered_species['unique_read_count'].max():,}")

        # Show top 10 species
        top_10 = filtered_species.nlargest(10, 'unique_read_count')
        print(f"\nTop 10 Species:")
        print(f"{'-'*40}")
        for idx, row in top_10.iterrows():
            print(f"{row['name']}: {row['unique_read_count']:,} unique reads")

        # Save detailed results to CSV
        results_df = filtered_species[['name', 'unique_read_count', 'read_count', 'percentage_abundance']].copy()
        results_df = results_df.sort_values('unique_read_count', ascending=False)
        results_df.columns = ['Species', 'Unique_Reads', 'Total_Reads', 'Percentage_Abundance']
        results_df.to_csv('classified_species_detailed.csv', index=False)
        print(f"\nDetailed results saved to: classified_species_detailed.csv")

        print(f"\n{'='*60}")
        print("INTERACTIVE PLOTS CREATED SUCCESSFULLY!")
        print(f"{'='*60}")
        print("1. classified_species_thin_bars.html - Standard thin bars")
        print("2. classified_species_very_thin_bars.html - Extra thin bars")
        print("3. classified_species_thin_bars_log.html - Thin bars with log scale")
        print("\nAll plots feature reduced bar width for better compactness.")

    except ImportError:
        print("Plotly is not installed. Installing it now...")
        import subprocess
        import sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "plotly"])
        print("Please run the script again now that Plotly is installed.")