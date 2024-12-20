#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""


from argparse import RawTextHelpFormatter
import logging
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string
import os
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from rich.progress import track
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def parser_harvest(parser):
    subparser = parser.add_parser(
        "harvest",
        description="""
Harvest predicted barcode start-end positions from barcode alignment count CSV files generated by
scarecrow seed.

Example:

scarecrow harveset BC1.csv BC2.csv BC3.csv --barcode_count 3 --min_distance 10
---
""",
        help="Harvest barcode start-end positions",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("barcodes", nargs="+", help="List of scarecrow barcode CSV output files")
    subparser.add_argument(
        "-o", "--out",
        metavar="barcodes_positions",
        help=("Path to output barcode positions file"),
        type=str,
        default="./barcode_positions.csv",
    )
    subparser.add_argument(
        "-b", "--barcode_count",
        metavar="barcode_count",
        help=("Number of barcodes expected [3]"),
        type=int,
        default=3,
    )    
    subparser.add_argument(
        "-d", "--min_distance",
        metavar="min_distance",
        help=("Minimum distance between peak start and end positions [10]"),
        type=int,
        default=10,
    )
    return subparser

def validate_harvest_args(parser, args):
    run_harvest(barcodes = args.barcodes,
                output_file = args.out,
                num_barcodes = args.barcode_count,
                min_distance = args.min_distance)

def run_harvest(barcodes, output_file, num_barcodes, min_distance):
    """
    Extract barcode positions from distribution of barcode alignments 
    Expects as an input a series of CSV files generated by scarecrow barcodes
    """
    
    # Global logger setup
    rndstring = generate_random_string()
    logfile = '{}_{}.{}'.format('./scarecrow_harvest', rndstring, 'log')
    logger = setup_logger(logfile)

    # Read barcode CSV files into a dataframe using pandas
    print(f"\033[34m\nImporting barcode data\033[0m")
    barcode_data = pd.concat([pd.read_csv(file, sep='\t') for file in barcodes], ignore_index=True)
    
    # Identify peaks
    results = get_barcode_peaks(barcode_data)
    top_peaks = top_peak_positions(results, num_barcodes, min_distance)
    print("\033[34m\nTop peaks\033[0m")
    print(f"{top_peaks}\n")

    # Plot
    pngfile = '{}_{}.{}'.format('./scarecrow_harvest', rndstring, 'png')
    print(f"\033[34m\nOutputting histogram to: \033[32m{pngfile}\033[0m")
    plot_peaks_faceted(barcode_data, outfile = pngfile)

    # Open file for writing output
    if output_file:
        top_peaks.to_csv(output_file, index = False)

    return 



def find_peaks_with_details(start_positions, end_positions, names):
    """
    Find peaks with detailed information about unique reads and names
    
    Args:
    - start_positions: Array of start positions
    - end_positions: Array of corresponding end positions
    - names: Array of corresponding read names
    
    Returns:
    - List of peak details: (start_position, end_position, unique_read_count, unique_name_fraction)
    """
    # Count frequency of each start position
    position_counts = pd.Series(start_positions).value_counts().sort_index()
    
    # Identify peaks in the frequency data
    scipy_peaks, _ = find_peaks(position_counts.values)
    
    # Extract peak positions with detailed information
    peaks_with_details = []
    for p in scipy_peaks:
        start = position_counts.index[p]
 
        # Filter to positions matching this start position
        mask = start_positions == start
        peak_names = names[mask]
        end = end_positions[mask][0]  # Get the corresponding end position
                
        # Calculate unique read count and name fraction
        peak_unique_names = np.unique(peak_names)
        read_count = len(peak_unique_names)
        
        # Ensure we don't divide by zero
        unique_names_count = len(np.unique(names))
        read_fraction = round(read_count / max(unique_names_count, 1),2)
              
        peaks_with_details.append((start, end, read_count, read_fraction))
    
    # Sort peaks by read_count in descending order
    peaks_with_details.sort(key=lambda x: x[2], reverse=True)
    
    return peaks_with_details


@log_errors
def get_barcode_peaks(barcode_data):
    """
    Read barcode CSV data in and identify peaks
    
    Args:
    - barcodes: List of CSV file paths
    """
    logger = logging.getLogger('scarecrow')
    
    # Identify barcode alignment peaks
    results = []
    barcode_groups = barcode_data.groupby(["read", "barcode_whitelist", "orientation"])
    #with tqdm(total = len(barcode_groups), desc = "Processing peaks") as pbar:
    for (read, barcode_whitelist, orientation), group in track(barcode_groups, total = len(barcode_groups)):
        # Find peaks with detailed information
        peaks = find_peaks_with_details(
            group["start"].values, 
            group["end"].values, 
            group["name"].values
        )
        
        # Add some debug logging
        logger.info(f"Peaks for {read} {barcode_whitelist} {orientation}:")
        for peak in peaks:
            logger.info(f"  Peak: start={peak[0]}, end={peak[1]}, read_count={peak[2]}, read_fraction={peak[3]}")
        
        # Store the results
        results.append({
            "read": read,
            "barcode_whitelist": barcode_whitelist,
            "orientation": orientation,
            "peaks": [
                {
                    "start": peak[0], 
                    "end": peak[1], 
                    "read_count": peak[2],
                    "read_fraction": peak[3]
                } for peak in peaks
            ]
        })
#        pbar.update(1)

    # Sort the results by the highest unique read count
    sorted_results = sorted(
        results,
        key=lambda x: max([peak["read_count"] for peak in x["peaks"]] or [0]),
        reverse=True
    )

    return sorted_results


def top_peak_positions(results, num_barcodes, min_distance):
    """
    Select top peaks, ensuring non-overlapping and minimum distance
    
    Args:
    - results: List of barcode peak results
    - num_barcodes: Number of barcodes to select
    - min_distance: Minimum distance between peaks
    """
    # Flatten the results to include individual peaks
    flattened_results = []
    for result in results:
        for peak in result["peaks"]:
            flattened_results.append({
                "start": peak["start"],
                "end": peak["end"],
                "read_count": peak["read_count"],
                "read_fraction": peak["read_fraction"],
                "barcode_whitelist": result["barcode_whitelist"],
                "read": result["read"],
                "orientation": result["orientation"]
            })

    # Convert to a DataFrame
    flattened_df = pd.DataFrame(flattened_results)

    # Group by position and barcode_whitelist, then aggregate
    grouped = flattened_df.groupby(
        ["read", "start", "end", "orientation", "barcode_whitelist"], 
        as_index=False
    ).agg({
        "read_count": "sum",
        "read_fraction": "mean"
    })

    # Sort by unique read count (descending) and position (ascending)
    sorted_grouped = grouped.sort_values(by=["read_count", "start"], ascending=[False, True])

    # Select top peaks with non-overlapping and minimum distance constraints
    top_positions = []
    for _, row in sorted_grouped.iterrows():
        # Check if this peak is far enough from all previously selected peaks
        if not any(
            abs(row["start"] - selected_peak["end"]) < min_distance or (row["start"] == selected_peak["start"])
            for selected_peak in top_positions
        ):
            top_positions.append(row)
            
            # Stop if we've found the desired number of peaks
            if len(top_positions) == num_barcodes:
                break

    # Create a DataFrame of the results
    top_positions_df = pd.DataFrame(top_positions).sort_values(by="start", ascending=True)

    return top_positions_df

def plot_peaks_faceted(barcode_data, save_fig=True, outfile='plot_faceted.png', dpi:int = 300):
    """
    Function to generate histogram of barcode start positions along each read in each barcode orientation
    """
   
    def plot_facet(data):
        g = sns.FacetGrid(data, col="read", row="barcode_whitelist", hue="read", margin_titles=True, height=3, aspect=3)
        g.map(sns.histplot, "start", binwidth=1, kde=False)
        g.set_axis_labels("Barcode start position", "Count", fontsize = 20)
        g.set_titles(col_template = "{col_name}", row_template = "{row_name}", size = 20)
        g.set(xlim=(barcode_data['start'].min(), barcode_data['start'].max()))
        for ax in g.axes.flat:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            if ax.texts:
                txt = ax.texts[0]
                ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1],
                        txt.get_text(), transform=ax.transAxes, va='center',
                        fontsize = 20)
                ax.texts[0].remove()
        return g
    
    # Create the 'forward' plot
    g1 = plot_facet(barcode_data[barcode_data['orientation'] == 'forward'])
    # Save the 'forward' plot temporarily
    g1.savefig('{}.{}'.format(outfile, 'f.png'), dpi = dpi)

    # Create the 'reverse' plot
    g2 = plot_facet(barcode_data[barcode_data['orientation'] == 'reverse']) 
    # Save the 'reverse' plot temporarily
    g2.savefig('{}.{}'.format(outfile, 'r.png'), dpi = dpi)

    # Combine the two plots side by side using matplotlib
    fig, axes = plt.subplots(1, 2, figsize=(11, 7))

    # Load the saved plots and display them
    forward_img = plt.imread('{}.{}'.format(outfile, 'f.png'))
    os.remove('{}.{}'.format(outfile, 'f.png'))
    reverse_img = plt.imread('{}.{}'.format(outfile, 'r.png'))
    os.remove('{}.{}'.format(outfile, 'r.png'))
    
    axes[0].imshow(forward_img)
    axes[0].axis('off')
    axes[0].set_title("Forward", fontsize = 6)

    axes[1].imshow(reverse_img)
    axes[1].axis('off')
    axes[1].set_title("Reverse", fontsize = 6)

    # Finalize and save the combined plot
    plt.savefig(outfile, dpi = dpi)
        
    return None