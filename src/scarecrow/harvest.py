# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

from argparse import RawTextHelpFormatter
import logging
import os
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.ticker import MaxNLocator
from functools import lru_cache
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


class ConservedRegion:
    """Represents a conserved region in a read sequence"""

    def __init__(
        self,
        file_index: int,  # Updated: Use file_index instead of read
        start: int,
        end: int,
        sequence: str,
        median_frequency: float,
    ):
        self.file_index = file_index  # Updated: Use file_index instead of read
        self.start = start
        self.end = end
        self.sequence = sequence
        self.median_frequency = median_frequency

    def __repr__(self) -> str:
        return f"ConservedRegion(file_index={self.file_index}, start={self.start}, end={self.end})"


class PeakAnalyzer:
    """Efficient peak analysis with caching and optimized computations"""

    def __init__(
        self,
        min_distance: int = 10,
        conserved_regions: Optional[List[ConservedRegion]] = None,
    ):
        self.min_distance = min_distance
        self._peak_cache = {}
        self.conserved_regions = conserved_regions or []

    @classmethod
    def from_conserved_file(
        cls, min_distance: int = 10, conserved_file: Optional[str] = None
    ) -> "PeakAnalyzer":
        """Create PeakAnalyzer instance with conserved regions loaded from file"""
        conserved_regions = []
        if conserved_file:
            df = pd.read_csv(conserved_file, sep="\t")
            conserved_regions = [
                ConservedRegion(
                    file_index=row["file_index"],  # Updated: Use file_index
                    start=row["start"],
                    end=row["end"],
                    sequence=row["sequence"],
                    median_frequency=row["median_frequency"],
                )
                for _, row in df.iterrows()
            ]
        return cls(min_distance=min_distance, conserved_regions=conserved_regions)

    @lru_cache(maxsize=128)
    def _compute_position_counts(self, positions: Tuple[int]) -> pd.Series:
        """Cached computation of position counts"""
        return pd.Series(positions).value_counts().sort_index()

    def _overlaps_conserved_region(self, file_index: int, start: int, end: int) -> bool:
        """Check if a peak overlaps with any conserved region"""
        for region in self.conserved_regions:
            if region.file_index == file_index and not (end < region.start or start > region.end):
                return True
        return False

    def find_peaks_optimized(
        self,
        start_positions: np.ndarray,
        end_positions: np.ndarray,
        names: np.ndarray,
        barcodes: np.ndarray,
        file_index: int,
    ) -> List[Tuple]:
        """Optimized peak finding with support for single peaks."""
        # Create cache key
        cache_key = hash((tuple(start_positions), tuple(end_positions), file_index))
        if cache_key in self._peak_cache:
            return self._peak_cache[cache_key]

        # Calculate position counts using cached method
        position_counts = self._compute_position_counts(tuple(start_positions))

        # Handle case where there is only one unique start position
        if len(position_counts) == 1:
            start = position_counts.index[0]
            end = end_positions[start_positions == start][0]
            count = position_counts.values[0]

            # Calculate barcode diversity
            unique_barcodes = len(np.unique(barcodes[start_positions == start]))
            total_barcodes = len(barcodes[start_positions == start])
            barcode_diversity = round(unique_barcodes / total_barcodes, 4)

            # Return the single peak
            peaks_with_details = [
                (start, end, count, 1.0, barcode_diversity)  # read_fraction = 1.0 for single peak
            ]
            self._peak_cache[cache_key] = peaks_with_details
            return peaks_with_details

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
            if self._overlaps_conserved_region(file_index, start, end):
                continue

            # Vectorized operations for peak statistics
            peak_unique_names = len(np.unique(peak_names))
            read_fraction = round(peak_unique_names / max(unique_names_count, 1), 2)

            # Calculate barcode diversity
            unique_barcodes = len(np.unique(peak_barcodes))
            total_barcodes = len(peak_barcodes)
            barcode_diversity = round(unique_barcodes / total_barcodes, 4)

            # Append results
            peaks_with_details.append(
                (start, end, peak_unique_names, read_fraction, barcode_diversity)
            )

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
    def _create_facet_plot(
        data: pd.DataFrame, ylim: Tuple[float, float]
    ) -> sns.FacetGrid:
        """Optimized facet plot creation"""
        g = sns.FacetGrid(
            data,
            col="read",
            row="barcode_whitelist",
            hue="read",
            margin_titles=True,
            height=3,
            aspect=3,
        )

        # Optimize plot rendering
        with plt.style.context("fast"):
            g.map(sns.histplot, "start", binwidth=1, kde=False)

        g.set_axis_labels("Barcode start position", "Count", fontsize=20)
        g.set_titles(col_template="{col_name}", row_template="{row_name}", size=20)
        g.set(xlim=(data["start"].min(), data["start"].max()))

        if ylim is not None:
            g.set(ylim=ylim)

        # Optimize axis rendering
        for ax in g.axes.flat:
            ax.tick_params(axis="x", labelsize=16)
            ax.tick_params(axis="y", labelsize=16)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
            if ax.texts:
                txt = ax.texts[0]
                ax.text(
                    txt.get_unitless_position()[0],
                    txt.get_unitless_position()[1],
                    txt.get_text(),
                    transform=ax.transAxes,
                    va="center",
                    fontsize=20,
                )
                ax.texts[0].remove()

        return g


class HarvestOptimizer:
    """Main class for optimized barcode harvesting"""

    def __init__(self, min_distance: int = 10, conserved_file: Optional[str] = None):
        self.peak_analyzer = PeakAnalyzer.from_conserved_file(
            min_distance, conserved_file
        )
        self.plotter = OptimizedPlotter()
        self.logger = logging.getLogger("scarecrow")

    def _overlaps_existing_peak(
        self,
        start: int,
        file_index: int,
        orientation: str,
        selected_peaks: List[Dict],
        min_distance: int,
    ) -> bool:
        """Check if a peak overlaps with any already selected peak across all whitelists"""
        for peak in selected_peaks:
            if (
                peak["file_index"] == file_index
                and peak["orientation"] == orientation
                and abs(start - peak["start"]) < min_distance
            ):
                return True
        return False

    def _select_peaks_for_whitelist(
        self,
        peaks_df: pd.DataFrame,
        num_barcodes: int,
        min_distance: int,
        existing_peaks: List[Dict],
    ) -> pd.DataFrame:
        """Select top peaks for a specific whitelist while avoiding overlap with existing peaks"""
        sorted_df = peaks_df.sort_values(
            ["read_count", "start"], ascending=[False, True]
        )
        selected_peaks = []

        for _, row in sorted_df.iterrows():
            # Check if peak overlaps with any existing peak across all whitelists
            if not self._overlaps_existing_peak(
                row["start"],
                row["file_index"],
                row["orientation"],
                existing_peaks,
                min_distance,
            ):
                selected_peaks.append(row.to_dict())
                existing_peaks.append(
                    row.to_dict()
                )  # Add to existing peaks to prevent future overlap
                if len(selected_peaks) == num_barcodes:
                    break

        return pd.DataFrame(selected_peaks)

    def _select_top_peaks(
        self, df: pd.DataFrame, num_barcodes: int, min_distance: int
    ) -> pd.DataFrame:
        """
        Select the best peak globally for each barcode_whitelist across all file_index values.
        Ensure that each peak is only assigned to one barcode_whitelist.
        """
        # Flatten peaks into a dataframe
        peak_data = []
        for _, row in df.iterrows():
            for peak in row["peaks"]:
                peak_data.append(
                    {
                        "barcode_whitelist": row["barcode_whitelist"],
                        "file_index": row["file_index"],
                        "file": row["file"],
                        "orientation": row["orientation"],
                        **peak,  # Unpack peak details (start, end, read_count, etc.)
                    }
                )

        if not peak_data:
            return pd.DataFrame()

        peaks_df = pd.DataFrame(peak_data)

        # Sort all peaks by read_count (descending) and start position (ascending)
        sorted_peaks = peaks_df.sort_values(
            ["read_count", "start"], ascending=[False, True]
        )

        # Track selected peaks and assigned whitelists
        selected_peaks = []
        assigned_whitelists = set()  # Track whitelists that have already been assigned
        assigned_file_indices = defaultdict(list)  # Track peaks assigned to each file_index

        for _, peak in sorted_peaks.iterrows():
            whitelist = peak["barcode_whitelist"]
            file_index = peak["file_index"]
            start = peak["start"]
            end = peak["end"]

            # Skip if the whitelist has already been assigned a peak
            if whitelist in assigned_whitelists:
                continue

            # Check if this peak overlaps with any already selected peak in the same file_index
            overlaps = False
            for existing_peak in assigned_file_indices[file_index]:
                if not (end < existing_peak["start"] or start > existing_peak["end"]):
                    overlaps = True
                    break

            if not overlaps:
                # Add this peak to the selected peaks
                selected_peaks.append(peak)
                assigned_whitelists.add(whitelist)  # Mark this whitelist as assigned
                assigned_file_indices[file_index].append(
                    {"start": start, "end": end}
                )

        # Convert selected peaks to a DataFrame
        if selected_peaks:
            result_df = pd.DataFrame(selected_peaks)
            return result_df
        return pd.DataFrame()

    def _verify_no_overlaps(self, df: pd.DataFrame, min_distance: int = 10) -> bool:
        """Verify that no peaks overlap within each file_index and orientation group."""
        # Group by file_index and orientation
        for (file_index, orientation), group in df.groupby(["file_index", "orientation"]):
            starts = group["start"].values
            ends = group["end"].values

            # Check for overlaps within this group
            for i in range(len(starts)):
                for j in range(i + 1, len(starts)):
                    # Check if the peaks overlap
                    if not (ends[i] < starts[j] or starts[i] > ends[j]):
                        # Peaks overlap
                        self.logger.error(
                            f"Found overlapping peaks in file {file_index}, orientation {orientation}: "
                            f"Peak 1 (start={starts[i]}, end={ends[i]}) overlaps with "
                            f"Peak 2 (start={starts[j]}, end={ends[j]})"
                        )
                        return False

        return True

    @log_errors
    def process_barcode_data(
        self, barcode_files: List[str], num_barcodes: int = 1, min_distance: int = 10
    ) -> pd.DataFrame:
        """Process barcode data and select unique top peaks per whitelist"""
        # Read data efficiently
        dfs = []
        for file in barcode_files:
            df = pd.read_csv(
                file,
                sep="\t",
                usecols=[
                    "file_index",
                    "file",
                    "read_name",
                    "seqlen",
                    "barcode_whitelist",
                    "barcode",
                    "orientation",
                    "start",
                    "end",
                ],
            )
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
        grouped = barcode_data.groupby(["file_index", "file", "barcode_whitelist", "orientation"])

        for (file_index, file, barcode_whitelist, orientation), group in grouped:

            # Convert to numpy arrays for faster processing
            start_arr = group["start"].values
            end_arr = group["end"].values
            names_arr = group["read_name"].values
            barcodes_arr = group["barcode"].values

            peaks = self.peak_analyzer.find_peaks_optimized(
                start_arr, end_arr, names_arr, barcodes_arr, file_index
            )

            results.append(
                {
                    "file_index": file_index,
                    "file": file,
                    "barcode_whitelist": barcode_whitelist,
                    "orientation": orientation,
                    "peaks": [
                        {
                            "start": p[0],
                            "end": p[1],
                            "read_count": p[2],
                            "read_fraction": p[3],
                        }
                        for p in peaks
                    ],
                }
            )

        # Sort efficiently using key function
        results.sort(
            key=lambda x: max((p["read_count"] for p in x["peaks"]), default=0),
            reverse=True,
        )
        return results


def parser_harvest(parser):
    subparser = parser.add_parser(
        "harvest",
        description="""
Generate barcode profiles by harvesting barcode matches identified by scarecrow seed.

Example:

scarecrow harveset BC1.csv BC2.csv BC3.csv \n\t--barcode_count 1 \n\t--min_distance 10 \n\t--conserved conserved.tsv \n\t--out barcode_positions.csv
---
""",
        help="Generate barcode profiles from barcode matches identified by scarecrow seed",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "barcodes", nargs="+", help="List of scarecrow barcode CSV output files"
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar="<file>",
        help=("Path to output barcode positions file"),
        type=str,
        default="./barcode_positions.csv",
    )
    subparser.add_argument(
        "-n",
        "--barcode_count",
        metavar="<int>",
        help=("Number of barcodes (peaks) to return per barcode:whitelist [1]"),
        type=int,
        default=1,
    )
    subparser.add_argument(
        "-d",
        "--min_distance",
        metavar="<int>",
        help=("Minimum distance between peak start and end positions [10]"),
        type=int,
        default=10,
    )
    subparser.add_argument(
        "-c",
        "--conserved",
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
    logfile = "{}_{}.{}".format(
        "./scarecrow_sam2fastq", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_harvest(
        barcodes=args.barcodes,
        output_file=args.out,
        num_barcodes=args.barcode_count,
        min_distance=args.min_distance,
        conserved_file=args.conserved,
    )


@log_errors
def run_harvest(
    barcodes: List[str],
    output_file: str,
    num_barcodes: int = 1,
    min_distance: int = 10,
    conserved_file: Optional[str] = None,
) -> None:
    """
    Optimized main function for harvesting barcode positions
    """
    logger = logging.getLogger("scarecrow")

    # Initialize optimizer
    optimizer = HarvestOptimizer(min_distance, conserved_file)

    # Process data
    results = optimizer.process_barcode_data(barcodes, num_barcodes, min_distance)

    # Log results by whitelist
    for whitelist, group in results.groupby("barcode_whitelist"):
        logger.info(f"\nPeaks for whitelist {whitelist} (n = {len(group)}):\n{group}")

    # Load conserved regions for plotting if file exists
    conserved_regions = None
    if conserved_file:
        logger.info(
            f"Excluded peaks overlapping with conserved regions from: {conserved_file}"
        )
        df = pd.read_csv(conserved_file, sep="\t")
        conserved_regions = [
            ConservedRegion(
                file_index=row["file_index"],
                start=row["start"],
                end=row["end"],
                sequence=row["sequence"],
                median_frequency=row["median_frequency"],
            )
            for _, row in df.iterrows()
        ]

    # Generate output
    if output_file:
        logger.info(f"Selected peaks written to: {output_file}")
        results.to_csv(output_file, index=False)

    # Generate plot with optimized settings
    logger.info(f"{results}")
    pngfile = Path(output_file).with_suffix('')
    plot_peaks_optimized(
        pd.concat([pd.read_csv(f, sep="\t") for f in barcodes], ignore_index=True),
        outfile=pngfile,
        conserved_regions=conserved_regions,
        selected_peaks=results,
    )
    logger.info("Finished!")
    return results


@log_errors
def plot_peaks_optimized(
    barcode_data: pd.DataFrame,
    outfile: str = "plot_faceted",
    dpi: int = 300,
    conserved_regions: Optional[List[ConservedRegion]] = None,
    selected_peaks: Optional[pd.DataFrame] = None,
) -> None:
    """Optimized plotting function with conserved region highlighting and peak annotations"""

    def _add_conserved_regions(ax, data, ylim):
        """Helper function to add conserved region highlighting"""
        if conserved_regions:
            file_index = data["file_index"].iloc[0] if not data.empty else None
            if file_index:
                for region in conserved_regions:
                    if region.file_index == file_index:  # Only add regions for matching read
                        # Add shaded region for conserved sequence
                        ax.axvspan(
                            region.start,
                            region.end,
                            alpha=0.2,
                            color="red",
                            label="Conserved Region",
                        )

    def _add_peak_annotations(ax, data, ylim):
        """Helper function to add peak annotations"""
        if selected_peaks is not None and not data.empty:
            read = data["file_index"].iloc[0]
            orientation = data["orientation"].iloc[0]
            whitelist = data["barcode_whitelist"].iloc[0]

            # Filter peaks for this specific read and orientation
            peaks = selected_peaks[
                (selected_peaks["file_index"] == read)
                & (selected_peaks["orientation"] == orientation)
            ]

            # Log filtered peaks for debugging
            # logger.info(f"Adding peak annotations for read={read}, orientation={orientation}, whitelist={whitelist}")
            # logger.info(f"Filtered peaks:\n{peaks}")

            for _, peak in peaks.iterrows():
                # Check if the peak's whitelist matches the current whitelist
                if peak["barcode_whitelist"] == whitelist:
                    # Add shaded region for peak
                    ax.axvspan(
                        peak["start"],
                        peak["end"],
                        alpha=0.2,
                        color="blue",
                        label="Peak",
                    )
                    # Log peak positions for debugging
                    # logger.info(f"Added peak annotation: start={peak['start']}, end={peak['end']}")
                # else:
                # logger.info(f"Skipping peak (whitelist mismatch): start={peak['start']}, end={peak['end']}, whitelist={peak['barcode_whitelist']}")

    def _create_enhanced_facet_plot(data, ylim):
        """Create facet plot with both conserved regions and peak annotations"""
        if data.empty:
            logger.warning("No data available for plotting. Skipping this group.")
            return None

        # Ensure the 'start' column exists and is valid
        if "start" not in data.columns or not pd.api.types.is_numeric_dtype(data["start"]):
            logger.warning("Invalid 'start' column. Skipping this group.")
            return None

        # Use `row` for barcode_whitelist to stack facets vertically
        g = sns.FacetGrid(
            data,
            row="barcode_whitelist",  # Facet by whitelist vertically
            hue="file_index",  # Use file_index for coloring
            margin_titles=True,
            height=3,
            aspect=3,
        )

        # Optimize plot rendering
        with plt.style.context("fast"):
            # Calculate the range of 'start' values
            start_range = data["start"].max() - data["start"].min()

            if start_range == 0:
                # If all 'start' values are the same, plot a single bar
                start_value = data["start"].iloc[0]
                count = len(data)

                # Create a bar plot instead of a histogram
                def plot_single_bar(data, **kwargs):
                    plt.bar(start_value, count, width = 0.8, align = "center")

                g.map(plot_single_bar, "start")
            else:
                # Use a default binwidth of 1, but ensure it's smaller than the range
                binwidth = min(0.5, start_range)
                bins = "auto"  # Let seaborn automatically determine the number of bins

                # Plot the histogram
                g.map(sns.histplot, "start", binwidth = binwidth, bins = bins, kde = False)

        g.set_axis_labels("Barcode start position", "Count", fontsize = 20)
        g.set_titles(row_template = "{row_name}", size = 20)  # Only row titles are needed
        g.fig.suptitle(f"File index: {file_index} ({file}), barcode orientation: {orientation}", y = 1.02)

        if ylim is not None:
            g.set(ylim = ylim)

        # Optimize axis rendering and add annotations
        for ax, (whitelist, subplot_data) in zip(
            g.axes.flat, data.groupby("barcode_whitelist")
        ):
            if not subplot_data.empty:
                # Determine per-subplot x-limits
                x_min, x_max = 1, subplot_data["start"].max()
                x_end = subplot_data.loc[
                    subplot_data["start"] == x_max, "seqlen"
                ].values[0]
                ax.set_xlim(x_min, x_end)  # Apply different x-limits per subplot
                _add_conserved_regions(ax, subplot_data, ylim)
                _add_peak_annotations(ax, subplot_data, ylim)

            ax.tick_params(axis="x", labelsize=16)
            ax.tick_params(axis="y", labelsize=16)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10))

            # Update subplot title if needed
            if ax.texts:
                txt = ax.texts[0]
                ax.text(
                    txt.get_unitless_position()[0],
                    txt.get_unitless_position()[1],
                    txt.get_text(),
                    transform=ax.transAxes,
                    va="center",
                    fontsize=20,
                )
                ax.texts[0].remove()

        return g

    # Get logger
    logger = logging.getLogger("scarecrow")

    # Calculate ylim
    max_count = (
        barcode_data.groupby(["file_index", "file", "orientation", "barcode_whitelist"])["start"]
        .value_counts()
        .max()
    )
    ylim = (0, max_count * 1.2)  # Add 20% padding for labels

    # Create plots using matplotlib's background renderer
    plt.switch_backend("Agg")

    # Separate data by file_index and orientation
    file_groups = barcode_data.groupby(["file_index", "file", "orientation"])

    # Create individual plots for each file and orientation
    logger.info("Generating plots of barcode count distributions")
    for (file_index, file, orientation), group in file_groups:
        if not group.empty:
            # Create plot for this file and orientation
            g = _create_enhanced_facet_plot(group, ylim)
            if g is not None:
                temp_file = f"{outfile}.file_index_{file_index}_{orientation}.png"
                g.savefig(temp_file, dpi=dpi, bbox_inches="tight")
                logger.info(f"- {temp_file}")
