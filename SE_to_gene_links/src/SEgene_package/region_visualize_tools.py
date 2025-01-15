from matplotlib.ticker import ScalarFormatter
from typing import Optional
import matplotlib.pyplot as plt
from pybedtools import BedTool


def plot_stacked_reads_bed(
    bed: BedTool,
    chrom: str,
    start: int,
    end: int,
    title: Optional[str] = None,
    xlabel: str = "Genomic Position",
    ylabel: str = "Read Stack",
    color: str = "blue",
    grid: bool = True,
    save_path: Optional[str] = None,
    save_region_bed: Optional[str] = None,
    save_full_bed: Optional[str] = None,
    show_plot: bool = True
) -> None:

    # Convert BED data to DataFrame
    bed_df = bed.to_dataframe()

    # Save full BED if requested
    if save_full_bed:
        bed.saveas(save_full_bed)
        print(f"Full BED saved as {save_full_bed}.")

    # Filter the data for the specified region
    region_data = bed_df[(bed_df['chrom'] == chrom) & 
                         (bed_df['start'] >= start) & 
                         (bed_df['end'] <= end)]

    if region_data.empty:
        print(f"No data found for the specified region: {chrom}:{start}-{end}")
        return

    # Save region BED if requested
    if save_region_bed:
        region_bed = BedTool.from_dataframe(region_data)
        region_bed.saveas(save_region_bed)
        print(f"Region BED saved as {save_region_bed}.")

    # Assign Y-axis positions to avoid overlaps
    region_data = region_data.sort_values(by=['start'])
    region_data['y'] = range(len(region_data))

    # Create the plot
    plt.figure(figsize=(15, 8))
    plt.hlines(region_data['y'], region_data['start'], region_data['end'], color=color)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plot_title = title if title else f"Stacked Read Plot for {chrom}:{start}-{end}"
    plt.title(plot_title, fontsize=18)
    plt.ylim(-1, len(region_data))
    
    # Configure axes
    plt.gca().xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    plt.ticklabel_format(style='plain', axis='x')

    # Add grid lines
    if grid:
        plt.grid(axis='x', linestyle='--', alpha=0.7)

    # Save the plot if a path is specified
    if save_path:
        plt.savefig(save_path, format='svg')
        print(f"Plot saved as {save_path}.")

    # Show the plot if requested
    if show_plot:
        plt.show()

    # Close the plot to free memory
    plt.clf()
    plt.close()