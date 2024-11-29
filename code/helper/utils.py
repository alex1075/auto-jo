from code.config.markers import markers
from fuzzywuzzy import fuzz

# Normalize channel names
def normalize_channel_name(channel):
    return channel.replace(" ", "").replace("*-A", "").replace("-A", "")


def prepare_marker_channels(sample, markers, output_file="channel_list.txt"):
    """
    Prepares a mapping of markers to channels by normalising channel names and performing fuzzy matching.

    Args:
        sample (fk.Sample): The FlowKit sample object.
        markers (dict): Dictionary mapping marker names to their expected channel names.
        output_file (str): File path to save the list of channels.

    Returns:
        dict: A dictionary mapping markers to their best-matched channels.
    """
    # Extract channel labels from the sample
    channels = sample.pnn_labels

    # Save channel list to a file
    with open(output_file, "w") as f:
        for channel in channels:
            f.write(channel + "\n")

    # Normalise channel names
    def normalize_channel_name(channel):
        return channel.replace(" ", "").replace("*-A", "").replace("-A", "")

    normalized_channels = {normalize_channel_name(channel): channel for channel in channels}

    # Match markers to channels using fuzzy matching
    marker_channels = {}
    for marker, marker_channel in markers.items():
        best_match = None
        best_score = 0
        normalized_marker_channel = normalize_channel_name(marker_channel)
        for normalized_channel, original_channel in normalized_channels.items():
            if normalized_channel == normalized_marker_channel:
                best_match = original_channel
                break
            score = fuzz.token_sort_ratio(normalized_channel, normalized_marker_channel)
            if score > best_score:
                best_score = score
                best_match = original_channel
        marker_channels[marker] = best_match

    return marker_channels


def return_marker_channels(marker, marker_channels):
    return marker_channels[marker]

def get_marker_index(marker, channels):
    for i, channel in enumerate(channels):
        if marker in channel:
            return i
    return None

