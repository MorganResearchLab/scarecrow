#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from argparse import RawTextHelpFormatter
import logging
import os
from typing import List, Dict, Tuple
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from functools import lru_cache
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


class PeakAnalyzer:
    """Efficient peak analysis with caching and optimized computations"""
    
    def __init__(self, min_distance: int = 10):
        self.min_distance = min_distance
        self._peak_cache = {}
    
    @lru_cache(maxsize=128)
    def _compute_position_counts(self, positions: Tuple[int]) -> pd.Series:
        """Cached computation of position counts"""
        return pd.Series(positions).value_counts().sort_index()
    
    def find_peaks_optimized(self, 
                           start_positions: np.ndarray, 
                           end_positions: np.ndarray, 
                           names: np.ndarray,
                           barcodes: np.ndarray) -> List[Tuple]:
        """Optimized peak finding with cached computations"""
        # Create cache key
        cache_key = hash((tuple(start_positions), tuple(end_positions)))
        if cache_key in self._peak_cache:
            return self._peak_cache[cache_key]

        # Calculate position counts using cached method
        position_counts = self._compute_position_counts(tuple(start_positions))
        
        # Optimize peak finding by pre-allocating arrays
        padded_values = np.zeros(len(position_counts) + 2)
        padded_values[1:-1] = position_counts.values
        
        # Find peaks efficiently
        peak_indices, _ = find_peaks(padded_values)
        peak_indices -= 1  # Adjust for padding
        
        # Process peaks efficiently using numpy operations
        peaks_with_details = []
        unique_names_count = len(np.unique(names))
        
        for p in peak_indices:
            start = position_counts.index[p]
            mask = start_positions == start
            peak_names = names[mask]
            peak_barcodes = barcodes[mask]
            end = end_positions[mask][0]
            
            # Vectorized operations for peak statistics
            peak_unique_names = len(np.unique(peak_names))
            read_fraction = round(peak_unique_names / max(unique_names_count, 1), 2)
            
            # Calculate barcode diversity
            unique_barcodes = len(np.unique(peak_barcodes))
            total_barcodes = len(peak_barcodes)
            barcode_diversity = round(unique_barcodes / total_barcodes, 4)

            # Append results
            peaks_with_details.append((start, end, peak_unique_names, read_fraction, barcode_diversity))
        
        # Sort peaks efficiently
        peaks_with_details.sort(key=lambda x: x[2], reverse=True)
        
        # Cache results
        self._peak_cache[cache_key] = peaks_with_details
        return peaks_with_details

class OptimizedPlotter:
    """Efficient plotting with caching and optimized rendering"""
    
    def __init__(self):
        self._plot_cache = {}
    
    @staticmethod
    def _create_facet_plot(data: pd.DataFrame, ylim: Tuple[float, float]) -> sns.FacetGrid:
        """Optimized facet plot creation"""
        g = sns.FacetGrid(
            data, 
            col="read", 
            row="barcode_whitelist", 
            hue="read",
            margin_titles=True, 
            height=3, 
            aspect=3
        )
        
        # Optimize plot rendering
        with plt.style.context('fast'):
            g.map(sns.histplot, "start", binwidth=1, kde=False)
        
        g.set_axis_labels("Barcode start position", "Count", fontsize=20)
        g.set_titles(col_template="{col_name}", row_template="{row_name}", size=20)
        g.set(xlim=(data['start'].min(), data['start'].max()))
        
        if ylim is not None:
            g.set(ylim=ylim)
        
        # Optimize axis rendering
        for ax in g.axes.flat:
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            if ax.texts:
                txt = ax.texts[0]
                ax.text(
                    txt.get_unitless_position()[0],
                    txt.get_unitless_position()[1],
                    txt.get_text(),
                    transform=ax.transAxes,
                    va='center',
                    fontsize=20
                )
                ax.texts[0].remove()
        
        return g

class HarvestOptimizer:
    """Main class for optimized barcode harvesting"""
    
    def __init__(self, min_distance: int = 10):
        self.peak_analyzer = PeakAnalyzer(min_distance)
        self.plotter = OptimizedPlotter()
        self.logger = logging.getLogger('scarecrow')
    
    @log_errors
    def process_barcode_data(self, 
                           barcode_files: List[str], 
                           num_barcodes: int, 
                           min_distance: int) -> pd.DataFrame:
        """Optimized barcode data processing"""
        # Read data efficiently
        dfs = []
        for file in barcode_files:
            df = pd.read_csv(file, sep='\t', usecols=['read', 'name', 'barcode_whitelist', 'barcode', 'orientation', 'start', 'end'])
            dfs.append(df)
        
        barcode_data = pd.concat(dfs, ignore_index=True)
        
        # Process peaks efficiently
        results = self._get_barcode_peaks(barcode_data)
        return self._select_top_peaks(results, num_barcodes, min_distance)
    
    def _get_barcode_peaks(self, barcode_data: pd.DataFrame) -> List[Dict]:
        """Optimized peak detection"""
        results = []
        
        # Group data efficiently
        grouped = barcode_data.groupby(['read', 'barcode_whitelist', 'orientation'])
        
        for (read, barcode_whitelist, orientation), group in grouped:
            # Convert to numpy arrays for faster processing
            start_arr = group['start'].values
            end_arr = group['end'].values
            names_arr = group['name'].values
            barcodes_arr = group['barcode'].values
            
            peaks = self.peak_analyzer.find_peaks_optimized(start_arr, end_arr, names_arr, barcodes_arr)
            
            results.append({
                'read': read,
                'barcode_whitelist': barcode_whitelist,
                'orientation': orientation,
                'peaks': [
                    {
                        'start': p[0],
                        'end': p[1],
                        'read_count': p[2],
                        'read_fraction': p[3],
                        'barcode_diversity': p[4]
                    } for p in peaks
                ]
            })
        
        # Sort efficiently using key function
        results.sort(key=lambda x: max((p['read_count'] for p in x['peaks']), default=0), reverse=True)
        return results
    
    def _select_top_peaks(self, 
                         results: List[Dict], 
                         num_barcodes: int, 
                         min_distance: int) -> pd.DataFrame:
        """Optimized peak selection"""
        # Flatten results efficiently
        flattened = [
            {
                'start': peak['start'],
                'end': peak['end'],
                'read_count': peak['read_count'],
                'read_fraction': peak['read_fraction'],
                'barcode_diversity': peak['barcode_diversity'],
                'barcode_whitelist': result['barcode_whitelist'],
                'read': result['read'],
                'orientation': result['orientation']
            }
            for result in results
            for peak in result['peaks']
        ]
        
        df = pd.DataFrame(flattened)
        
        # Efficient grouping and aggregation
        grouped = df.groupby(
            ['read', 'start', 'end', 'orientation', 'barcode_whitelist'],
            as_index=False
        ).agg({
            'barcode_diversity': 'median',
            'read_count': 'sum',
            'read_fraction': 'median'
        })
        
        sorted_df = grouped.sort_values(['read_count', 'start'], ascending=[False, True])
        
        # Efficient peak selection
        selected_peaks = []
        for _, row in sorted_df.iterrows():
            if not any(abs(row['start'] - peak['end']) < min_distance or 
                      (row['start'] == peak['start']) for peak in selected_peaks):
                selected_peaks.append(row)
                if len(selected_peaks) == num_barcodes:
                    break
        
        return pd.DataFrame(selected_peaks).sort_values('read_count', ascending=False)


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
    
@log_errors
def run_harvest(barcodes: List[str],
                output_file: str,
                num_barcodes: int = 3,
                min_distance: int = 10) -> None:
    """
    Optimized main function for harvesting barcode positions
    """
    # Setup logging
    logfile = f'./scarecrow_harvest_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    # Initialize optimizer
    optimizer = HarvestOptimizer(min_distance)
    
    # Process data
    results = optimizer.process_barcode_data(barcodes, num_barcodes, min_distance)
    logger.info(f"Peaks (n = {num_barcodes}) identified:\n{results}")

    # Generate output
    if output_file:
        results.to_csv(output_file, index=False)
    
    # Generate plot with optimized settings
    pngfile = f'./scarecrow_harvest_{generate_random_string()}.png'
    plot_peaks_optimized(pd.concat([pd.read_csv(f, sep='\t') for f in barcodes], ignore_index=True), 
                        outfile=pngfile)
    
    return results

@log_errors
def plot_peaks_optimized(barcode_data: pd.DataFrame, 
                        outfile: str = 'plot_faceted.png',
                        dpi: int = 300) -> None:
    """Optimized plotting function"""
    logger = logging.getLogger('scarecrow')

    plotter = OptimizedPlotter()
    
    # Calculate ylim efficiently
    max_count = barcode_data.groupby(['read', 'orientation', 'barcode_whitelist'])['start'].value_counts().max()
    ylim = (0, max_count)
    
    # Create plots efficiently using matplotlib's background renderer
    plt.switch_backend('Agg')
    
    # Generate plots for forward and reverse orientations
    forward_data = barcode_data[barcode_data['orientation'] == 'forward']
    reverse_data = barcode_data[barcode_data['orientation'] == 'reverse']
    
    g1 = plotter._create_facet_plot(forward_data, ylim)
    g2 = plotter._create_facet_plot(reverse_data, ylim)
    
    # Save temporary files with optimized compression
    temp_forward = f'{outfile}.f.png'
    temp_reverse = f'{outfile}.r.png'
    
    g1.savefig(temp_forward, dpi=dpi, bbox_inches='tight')
    g2.savefig(temp_reverse, dpi=dpi, bbox_inches='tight')
    
    # Create final combined plot
    fig, axes = plt.subplots(2, 1, figsize=(4, 4))
    
    # Load and display images efficiently
    for img_file, ax, title in [
        (temp_forward, axes[0], "Forward orientation"),
        (temp_reverse, axes[1], "Reverse orientation")
    ]:
        img = plt.imread(img_file)
        ax.imshow(img)
        ax.axis('off')
        ax.set_title(title, fontsize=4, loc='left')
        os.remove(img_file)
    
    # Save final plot efficiently
    fig.tight_layout(pad=0)
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"barcode count distribution: {outfile}")