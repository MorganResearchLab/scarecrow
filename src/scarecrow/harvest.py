#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from argparse import RawTextHelpFormatter
import logging
import os
from typing import List, Dict, Tuple, Optional
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from functools import lru_cache
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

class ConservedRegion:
    """Represents a conserved region in a read sequence"""
    def __init__(self, read: str, start: int, end: int, sequence: str, median_frequency: float):
        self.read = read
        self.start = start
        self.end = end
        self.sequence = sequence
        self.median_frequency = median_frequency

    def __repr__(self) -> str:
        return f"ConservedRegion(read={self.read}, start={self.start}, end={self.end})"

class PeakAnalyzer:
    """Efficient peak analysis with caching and optimized computations"""
    
    def __init__(self, min_distance: int = 10, conserved_regions: Optional[List[ConservedRegion]] = None):
        self.min_distance = min_distance
        self._peak_cache = {}
        self.conserved_regions = conserved_regions or []
    
    @classmethod
    def from_conserved_file(cls, min_distance: int = 10, conserved_file: Optional[str] = None) -> 'PeakAnalyzer':
        """Create PeakAnalyzer instance with conserved regions loaded from file"""
        conserved_regions = []
        if conserved_file:
            df = pd.read_csv(conserved_file, sep='\t')
            conserved_regions = [
                ConservedRegion(
                    read = row['read'],
                    start = row['start'],
                    end = row['end'],
                    sequence = row['sequence'],
                    median_frequency = row['median_frequency']
                )
                for _, row in df.iterrows()
            ]
        return cls(min_distance = min_distance, conserved_regions = conserved_regions)
    
    @lru_cache(maxsize=128)
    def _compute_position_counts(self, positions: Tuple[int]) -> pd.Series:
        """Cached computation of position counts"""
        return pd.Series(positions).value_counts().sort_index()
    
    def _overlaps_conserved_region(self, read: str, start: int, end: int) -> bool:
        """Check if a peak overlaps with any conserved region"""
        for region in self.conserved_regions:
            if (region.read == read and 
                not (end < region.start or start > region.end)):
                return True
        return False

    def find_peaks_optimized(self, 
                           start_positions: np.ndarray, 
                           end_positions: np.ndarray, 
                           names: np.ndarray,
                           barcodes: np.ndarray,
                           read: str) -> List[Tuple]:
        """Optimized peak finding with cached computations and conserved region filtering"""
        # Create cache key
        cache_key = hash((tuple(start_positions), tuple(end_positions), read))
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
            
            # Skip peaks that overlap with conserved regions
            if self._overlaps_conserved_region(read, start, end):
                continue
            
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
            g.map(sns.histplot, "start", binwidth = 1, kde=False)
        
        g.set_axis_labels("Barcode start position", "Count", fontsize = 20)
        g.set_titles(col_template="{col_name}", row_template="{row_name}", size = 20)
        g.set(xlim=(data['start'].min(), data['start'].max()))
        
        if ylim is not None:
            g.set(ylim = ylim)
        
        # Optimize axis rendering
        for ax in g.axes.flat:
            ax.tick_params(axis = 'x', labelsize = 16)
            ax.tick_params(axis = 'y', labelsize = 16)
            ax.xaxis.set_major_locator(MaxNLocator(nbins = 10))
            if ax.texts:
                txt = ax.texts[0]
                ax.text(
                    txt.get_unitless_position()[0],
                    txt.get_unitless_position()[1],
                    txt.get_text(),
                    transform = ax.transAxes,
                    va = 'center',
                    fontsize = 20
                )
                ax.texts[0].remove()
        
        return g

class HarvestOptimizer:
    """Main class for optimized barcode harvesting"""
    
    def __init__(self, min_distance: int = 10, conserved_file: Optional[str] = None):
        self.peak_analyzer = PeakAnalyzer.from_conserved_file(min_distance, conserved_file)
        self.plotter = OptimizedPlotter()
        self.logger = logging.getLogger('scarecrow')
    
    def _overlaps_existing_peak(self, 
                              start: int, 
                              read: str, 
                              orientation: str,
                              selected_peaks: List[Dict],
                              min_distance: int) -> bool:
        """Check if a peak overlaps with any already selected peak across all whitelists"""
        for peak in selected_peaks:
            if (peak['read'] == read and 
                peak['orientation'] == orientation and
                abs(start - peak['start']) < min_distance):
                return True
        return False
    
    def _select_peaks_for_whitelist(self, 
                                  peaks_df: pd.DataFrame, 
                                  num_barcodes: int, 
                                  min_distance: int,
                                  existing_peaks: List[Dict]) -> pd.DataFrame:
        """Select top peaks for a specific whitelist while avoiding overlap with existing peaks"""
        sorted_df = peaks_df.sort_values(['read_count', 'start'], ascending=[False, True])
        selected_peaks = []
        
        for _, row in sorted_df.iterrows():
            # Check if peak overlaps with any existing peak across all whitelists
            if not self._overlaps_existing_peak(row['start'], 
                                              row['read'], 
                                              row['orientation'],
                                              existing_peaks, 
                                              min_distance):
                selected_peaks.append(row.to_dict())
                existing_peaks.append(row.to_dict())  # Add to existing peaks to prevent future overlap
                if len(selected_peaks) == num_barcodes:
                    break
        
        return pd.DataFrame(selected_peaks)
    
    def _select_top_peaks(self, 
                         df: pd.DataFrame, 
                         num_barcodes: int, 
                         min_distance: int) -> pd.DataFrame:
        """Select top peaks for each whitelist separately, ensuring no peaks are shared"""
        # Flatten peaks into a dataframe
        peak_data = []
        for _, row in df.iterrows():
            for peak in row['peaks']:
                peak_data.append({
                    'barcode_whitelist': row['barcode_whitelist'],
                    'read': row['read'],
                    'orientation': row['orientation'],
                    **peak  # Unpack peak details
                })
        
        if not peak_data:
            return pd.DataFrame()
            
        peaks_df = pd.DataFrame(peak_data)
        
        # Sort whitelists by total read count to prioritize stronger signals
        whitelist_priorities = (peaks_df.groupby('barcode_whitelist')['read_count']
                              .sum()
                              .sort_values(ascending=False))
        
        # Select peaks for each whitelist in priority order
        selected_peaks = []
        all_selected_peaks = []  # Track all selected peaks across whitelists
        
        for whitelist in whitelist_priorities.index:
            whitelist_peaks = peaks_df[peaks_df['barcode_whitelist'] == whitelist]
            whitelist_selected = self._select_peaks_for_whitelist(
                whitelist_peaks, 
                num_barcodes, 
                min_distance,
                all_selected_peaks  # Pass all selected peaks to avoid overlap
            )
            selected_peaks.append(whitelist_selected)
        
        # Combine all selected peaks
        if selected_peaks:
            result_df = pd.concat(selected_peaks, ignore_index = True)
            # Double-check for any overlaps
            assert self._verify_no_overlaps(result_df, min_distance), "Found overlapping peaks"
            return result_df
        return pd.DataFrame()
    
    def _verify_no_overlaps(self, df: pd.DataFrame, min_distance: int = 10) -> bool:
        """Verify that no peaks overlap within each read and orientation"""
        for (read, orientation), group in df.groupby(['read', 'orientation']):
            starts = group['start'].values
            for i, start1 in enumerate(starts):
                for start2 in starts[i+1:]:
                    if abs(start1 - start2) < min_distance:
                        self.logger.error(f"Found overlapping peaks: {start1} and {start2} "
                                        f"in read {read}, orientation {orientation}")
                        return False
        return True
    
    @log_errors
    def process_barcode_data(self, 
                            barcode_files: List[str], 
                            num_barcodes: int = 1, 
                            min_distance: int = 10) -> pd.DataFrame:
        """Process barcode data and select unique top peaks per whitelist"""
        # Read data efficiently
        dfs = []
        for file in barcode_files:
            df = pd.read_csv(file, sep = '\t', 
                           usecols=['read', 'name', 'seqlen', 'barcode_whitelist', 
                                  'barcode', 'orientation', 'start', 'end'])
            dfs.append(df)
        
        barcode_data = pd.concat(dfs, ignore_index=True)
        
        # Process peaks efficiently
        peaks_by_group = self._get_barcode_peaks(barcode_data)
        peaks_df = pd.DataFrame(peaks_by_group)

        self.logger.info(f"Top 10 peaks\n{peaks_df.head(10)}")
        
        # Select unique top peaks for each whitelist
        return self._select_top_peaks(peaks_df, num_barcodes, min_distance)

    def _get_barcode_peaks(self, barcode_data: pd.DataFrame) -> List[Dict]:
        """Optimized peak detection with conserved region filtering"""
        results = []
        
        # Group data efficiently
        grouped = barcode_data.groupby(['read', 'barcode_whitelist', 'orientation'])
        
        for (read, barcode_whitelist, orientation), group in grouped:
            # Convert to numpy arrays for faster processing
            start_arr = group['start'].values
            end_arr = group['end'].values
            names_arr = group['name'].values
            barcodes_arr = group['barcode'].values
            
            peaks = self.peak_analyzer.find_peaks_optimized(
                start_arr, end_arr, names_arr, barcodes_arr, read
            )
            
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



def parser_harvest(parser):
    subparser = parser.add_parser(
        "harvest",
        description="""
Harvest predicted barcode start-end positions from barcode alignment count CSV files generated by
scarecrow seed.

Example:

scarecrow harveset BC1.csv BC2.csv BC3.csv \n\t--barcode_count 1 \n\t--min_distance 10 \n\t--conserved conserved.tsv \n\t--out barcode_positions.csv
---
""",
        help="Harvest barcode start-end positions",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("barcodes", nargs="+", help="List of scarecrow barcode CSV output files")
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to output barcode positions file"),
        type=str,
        default="./barcode_positions.csv",
    )
    subparser.add_argument(
        "-n", "--barcode_count",
        metavar="<int>",
        help=("Number of barcodes (peaks) to return per barcode:whitelist [1]"),
        type=int,
        default=1,
    )    
    subparser.add_argument(
        "-d", "--min_distance",
        metavar="<int>",
        help=("Minimum distance between peak start and end positions [10]"),
        type=int,
        default=10,
    )
    subparser.add_argument(
        "-c", "--conserved",
        metavar="<file>",
        help=("Path to TSV file containing conserved regions to exclude"),
        type=str,
        default=None,
    )
    return subparser

def validate_harvest_args(parser, args):
    """ 
    Validate arguments 
    """
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_sam2fastq', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_harvest(barcodes = args.barcodes,
                output_file = args.out,
                num_barcodes = args.barcode_count,
                min_distance = args.min_distance,
                conserved_file = args.conserved)
    
@log_errors
def run_harvest(barcodes: List[str],
                output_file: str,
                num_barcodes: int = 1,
                min_distance: int = 10,
                conserved_file: Optional[str] = None) -> None:
    """
    Optimized main function for harvesting barcode positions
    """
    logger = logging.getLogger('scarecrow')
    
    # Initialize optimizer
    optimizer = HarvestOptimizer(min_distance, conserved_file)
    
    # Process data
    results = optimizer.process_barcode_data(barcodes, num_barcodes, min_distance)
    
    # Log results by whitelist
    for whitelist, group in results.groupby('barcode_whitelist'):
        logger.info(f"\nPeaks for whitelist {whitelist} (n = {len(group)}):\n{group}")

    # Load conserved regions for plotting if file exists
    conserved_regions = None
    if conserved_file:
        logger.info(f"Excluded peaks overlapping with conserved regions from: {conserved_file}")
        df = pd.read_csv(conserved_file, sep = '\t')
        conserved_regions = [
            ConservedRegion(
                read = row['read'],
                start = row['start'],
                end = row['end'],
                sequence = row['sequence'],
                median_frequency = row['median_frequency']
            )
            for _, row in df.iterrows()
        ]

    # Generate output
    if output_file:
        logger.info(f"Selected peaks written to: {output_file}")
        results.to_csv(output_file, index = False)
    
    # Generate plot with optimized settings
    pngfile = f'./scarecrow_harvest_{generate_random_string()}.png'
    plot_peaks_optimized(
        pd.concat([pd.read_csv(f, sep = '\t') for f in barcodes], ignore_index = True), 
        outfile = pngfile,
        conserved_regions = conserved_regions,
        selected_peaks = results
    )
    
    return results

@log_errors
def plot_peaks_optimized(barcode_data: pd.DataFrame, 
                         outfile: str = 'plot_faceted.png',
                         dpi: int = 300,
                         conserved_regions: Optional[List[ConservedRegion]] = None,
                         selected_peaks: Optional[pd.DataFrame] = None) -> None:
    """Optimized plotting function with conserved region highlighting and peak annotations"""
    logger = logging.getLogger('scarecrow')

    plotter = OptimizedPlotter()
    
    # Calculate ylim
    max_count = barcode_data.groupby(['read', 'orientation', 'barcode_whitelist'])['start'].value_counts().max()
    ylim = (0, max_count * 1.2)  # Add 20% padding for labels
    
    # Create plots using matplotlib's background renderer
    plt.switch_backend('Agg')
    
    # Separate data into read1 and read2 for forward and reverse orientations
    read1_forward = barcode_data[(barcode_data['read'] == 'read1') & (barcode_data['orientation'] == 'forward')]
    read1_reverse = barcode_data[(barcode_data['read'] == 'read1') & (barcode_data['orientation'] == 'reverse')]
    read2_forward = barcode_data[(barcode_data['read'] == 'read2') & (barcode_data['orientation'] == 'forward')]
    read2_reverse = barcode_data[(barcode_data['read'] == 'read2') & (barcode_data['orientation'] == 'reverse')]

    def _add_conserved_regions(ax, data, ylim):
        """Helper function to add conserved region highlighting"""
        if conserved_regions:
            read = data['read'].iloc[0] if not data.empty else None
            if read:
                for region in conserved_regions:
                    if region.read == read:  # Only add regions for matching read
                        # Add shaded region for conserved sequence
                        ax.axvspan(region.start, region.end, 
                                   alpha=0.2, 
                                   color='red', 
                                   label='Conserved Region')

    def _add_peak_annotations(ax, data, ylim):
        """Helper function to add peak annotations"""
        if selected_peaks is not None and not data.empty:
            read = data['read'].iloc[0]
            orientation = data['orientation'].iloc[0]
            whitelist = data['barcode_whitelist'].iloc[0]
            
            # Filter peaks for this specific read and orientation
            peaks = selected_peaks[
                (selected_peaks['read'] == read) & 
                (selected_peaks['orientation'] == orientation)
            ]
            
            # Log filtered peaks for debugging
            #logger.info(f"Adding peak annotations for read={read}, orientation={orientation}, whitelist={whitelist}")
            #logger.info(f"Filtered peaks:\n{peaks}")
            
            for _, peak in peaks.iterrows():
                # Check if the peak's whitelist matches the current whitelist
                if peak['barcode_whitelist'] == whitelist:
                    # Add shaded region for peak
                    ax.axvspan(peak['start'], peak['end'], 
                            alpha=0.2, 
                            color='blue', 
                            label="Peak")
                    # Log peak positions for debugging
                    #logger.info(f"Added peak annotation: start={peak['start']}, end={peak['end']}")
                #else:
                    #logger.info(f"Skipping peak (whitelist mismatch): start={peak['start']}, end={peak['end']}, whitelist={peak['barcode_whitelist']}")
    
    def _create_enhanced_facet_plot(data, ylim):
        """Create facet plot with both conserved regions and peak annotations"""
        if data.empty:
            return None
        
        # Use `row` for barcode_whitelist to stack facets vertically
        g = sns.FacetGrid(
            data, 
            row="barcode_whitelist",  # Facet by whitelist vertically
            hue="read",
            margin_titles=True, 
            height=3, 
            aspect=3
        )
        
        # Optimize plot rendering
        with plt.style.context('fast'):
            g.map(sns.histplot, "start", binwidth=1, kde=False)
        
        g.set_axis_labels("Barcode start position", "Count", fontsize=20)
        g.set_titles(row_template="{row_name}", size=20)  # Only row titles are needed
        
        if ylim is not None:
            g.set(ylim=ylim)        

        # Optimize axis rendering and add annotations
        for ax, (whitelist, subplot_data) in zip(g.axes.flat, data.groupby('barcode_whitelist')):
            if not subplot_data.empty:
                # Determine per-subplot x-limits    
                x_min, x_max = subplot_data['start'].min(), subplot_data['start'].max()
                x_end = subplot_data.loc[subplot_data['start'] == x_max, 'seqlen'].values[0]                
                ax.set_xlim(x_min, x_end)  # Apply different x-limits per subplot
                _add_conserved_regions(ax, subplot_data, ylim)
                _add_peak_annotations(ax, subplot_data, ylim)
            
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            
            # Update subplot title if needed
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

    # Create individual plots for each read and orientation
    g1 = _create_enhanced_facet_plot(read1_forward, ylim)
    g2 = _create_enhanced_facet_plot(read1_reverse, ylim)
    g3 = _create_enhanced_facet_plot(read2_forward, ylim)
    g4 = _create_enhanced_facet_plot(read2_reverse, ylim)

    # Save temporary files with optimized compression
    temp_files = []
    for i, g in enumerate([g1, g2, g3, g4]):
        if g is not None:
            temp_file = f'{outfile}.temp_{i}.png'
            g.savefig(temp_file, dpi=dpi, bbox_inches='tight')
            temp_files.append(temp_file)

    # Create final combined plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # Load and display images in the correct order
    for i, (img_file, ax, title) in enumerate([
        (temp_files[0], axes[0, 0], "Read 1: forward"),
        (temp_files[1], axes[0, 1], "Read 1: reverse"),
        (temp_files[2], axes[1, 0], "Read 2: forward"),
        (temp_files[3], axes[1, 1], "Read 2: reverse")
    ]):
        if os.path.exists(img_file):
            img = plt.imread(img_file)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(title, fontsize=12, loc='left')
            os.remove(img_file)
    
    # Save final plot efficiently
    fig.tight_layout(pad=2.0)
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"Barcode count distribution with annotations: {outfile}")