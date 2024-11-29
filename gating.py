import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def gate_autofluorescence_vs_ssca(sample, autofluorescence_channel, ssca_channel, source="xform", show_hist=False):
    """
    Gate populations based on autofluorescence against SSC-A.

    Args:
        sample (fk.Sample): The FlowKit sample object.
        autofluorescence_channel (str): Channel name for autofluorescence.
        ssca_channel (str): Channel name for SSC-A.
        source (str): Source of the data ('orig', 'raw', 'comp', or 'xform'). Default is 'xform'.
        show_hist (bool): Whether to show histogram plots for debugging.

    Returns:
        dict: Dictionary containing gated populations as NumPy arrays.
    """
    # Map marker names to channel indices
    auto_index = sample.pnn_labels.index(autofluorescence_channel)
    ssca_index = sample.pnn_labels.index(ssca_channel)

    # Extract event data for the two channels
    data = sample.get_events(source=source)
    autofluorescence_data = data[:, auto_index]
    ssca_data = data[:, ssca_index]

    # Preprocess data: Remove invalid or non-positive values
    autofluorescence_data = np.where(autofluorescence_data > 0, autofluorescence_data, 1e-5)
    ssca_data = np.where(ssca_data > 0, ssca_data, 1e-5)

    # Histogram-based thresholding for autofluorescence
    hist, bin_edges = np.histogram(autofluorescence_data, bins=100)
    peaks, _ = find_peaks(hist, height=0.1 * max(hist), distance=5)

    # Determine thresholds for gating
    if len(peaks) >= 2:
        # Find valleys between peaks (local minima)
        valleys, _ = find_peaks(-hist, height=-0.1 * max(hist), distance=5)
        autofluorescence_threshold = bin_edges[valleys[0]]  # First valley as the threshold
    else:
        # Fallback to median if peaks not found
        autofluorescence_threshold = np.median(autofluorescence_data)

    # Debug: Show histogram
    if show_hist:
        plt.figure(figsize=(8, 4))
        plt.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor="k", alpha=0.7)
        plt.scatter(bin_edges[peaks], hist[peaks], color="red", label="Peaks")
        plt.axvline(autofluorescence_threshold, color="blue", linestyle="--", label="Threshold")
        plt.title("Autofluorescence Histogram")
        plt.legend()
        plt.show()

    # Gating logic
    high_auto = data[autofluorescence_data > autofluorescence_threshold]
    low_auto = data[autofluorescence_data <= autofluorescence_threshold]

    # Return gated populations
    populations = {
        "High Autofluorescence": high_auto,
        "Low Autofluorescence": low_auto,
    }

    return populations


def compute_threshold_histogram(data, bins=100, max_transform_attempts=5, show_hist=False):
    """
    Compute thresholds based on histogram peaks with axis transformations.

    Args:
        data (np.array): 1D array of marker intensities.
        bins (int): Number of bins for the histogram.
        max_transform_attempts (int): Maximum number of transformations to try.
        show_hist (bool): Whether to display histograms for debugging.

    Returns:
        list: List of thresholds (valleys between peaks) or fallback (median) if no peaks are found.
    """
    def transform_axis(data, attempt):
        """Apply a transformation to the data based on the attempt."""
        if attempt == 0:
            return data  # No transformation
        elif attempt == 1:
            return np.log1p(data)  # Logarithmic
        elif attempt == 2:
            return np.sqrt(data)  # Square root
        elif attempt == 3:
            return np.cbrt(data)  # Cube root
        elif attempt == 4:
            return np.expm1(data / max(data))  # Exponential (scaled)
        else:
            raise ValueError("Max transformation attempts exceeded.")

    # Pre-process data: remove NaN/Inf and replace non-positive values
    data = data[np.isfinite(data)]  # Remove NaN and Inf
    data = np.where(data > 0, data, 1e-5)  # Replace non-positive values with a small positive number

    for attempt in range(max_transform_attempts):
        # Transform the data
        transformed_data = transform_axis(data, attempt)

        # Compute histogram
        hist, bin_edges = np.histogram(transformed_data, bins=bins)

        # Find peaks (maxima) in the histogram
        peaks, _ = find_peaks(hist, height=0.1 * max(hist), distance=5)

        # If at least two peaks are found, stop transforming
        if len(peaks) >= 2:
            # Identify valleys (local minima) between peaks
            valleys, _ = find_peaks(-hist, height=-0.1 * max(hist), distance=5)

            # Filter valleys to include only those between the peaks
            thresholds = [
                bin_edges[valleys[i]] for i in range(len(valleys))
                if valleys[i] > peaks[0] and valleys[i] < peaks[-1]
            ]

            # Debugging: Plot the histogram if requested
            if show_hist:
                plt.figure(figsize=(8, 4))
                plt.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor="k", alpha=0.7)
                plt.scatter(bin_edges[peaks], hist[peaks], color="red", label="Peaks")
                plt.scatter(bin_edges[valleys], hist[valleys], color="blue", label="Valleys")
                plt.title(f"Histogram with {len(peaks)} Peaks (Attempt {attempt})")
                plt.legend()
                plt.show()

            return thresholds

    # Fallback: Use the median of the data as a threshold
    print("No clear peaks found after transformations. Falling back to median threshold.")
    return [np.median(data)]

def get_marker_populations(sample, marker1, marker2, source="xform"):
    """
    Function to extract populations based on two markers from a FlowKit sample.

    Args:
        sample (fk.Sample): A FlowKit sample object.
        marker1 (str): The name of the first marker.
        marker2 (str): The name of the second marker.
        source (str): Source of the data ('orig', 'raw', 'comp', or 'xform'). Default is 'xform'.

    Returns:
        dict: A dictionary containing the populations as NumPy arrays, labeled with f-strings.
    """
    # Map marker names to channel names
    channel1 = return_marker_channels(marker1)
    channel2 = return_marker_channels(marker2)

    # Debug: Ensure channels are correctly mapped
    print(f"Using channels: {channel1} for {marker1}, {channel2} for {marker2}")

    # Extract data for the two channels
    data = sample.get_events(source=source)

    # Find indices of the channels in the sample
    channel1_index = sample.pnn_labels.index(channel1)
    channel2_index = sample.pnn_labels.index(channel2)

    # Extract channel data
    channel1_data = data[:, channel1_index]
    channel2_data = data[:, channel2_index]

    # Compute thresholds
    marker1_thresholds = compute_threshold_histogram(channel1_data)
    marker2_thresholds = compute_threshold_histogram(channel2_data)

    # Use the first threshold (fallback logic for single threshold)
    marker1_threshold = marker1_thresholds[0] if marker1_thresholds else np.median(channel1_data)
    marker2_threshold = marker2_thresholds[0] if marker2_thresholds else np.median(channel2_data)

    # Create masks for populations
    populations = {
        f"{marker1}+ {marker2}-": data[(channel1_data > marker1_threshold) & (channel2_data <= marker2_threshold)],
        f"{marker1}+ {marker2}+": data[(channel1_data > marker1_threshold) & (channel2_data > marker2_threshold)],
        f"{marker1}- {marker2}+": data[(channel1_data <= marker1_threshold) & (channel2_data > marker2_threshold)],
        f"{marker1}- {marker2}-": data[(channel1_data <= marker1_threshold) & (channel2_data <= marker2_threshold)],
    }

    return populations