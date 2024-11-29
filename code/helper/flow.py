import flowkit as fk

def read_sample(path):
    """
    Function to read a FlowKit sample from a file.

    Args:
        path (str): The path to the FlowKit sample file.

    Returns:
        fk.Sample: The FlowKit sample object.
    """
    return fk.Sample(path)

def make_channel_list(sample, output_path):
    """
    Function to write the list of channels in a FlowKit sample to a file.

    Args:
        sample (fk.Sample): The FlowKit sample object.
        output_path (str): The path to the output file.
    """
    channels = sample.pnn_labels
    return channels

def normalize_channel_name(channel):
    return channel.replace(" ", "").replace("*-A", "").replace("-A", "")

def normalize_channels(channels):
    return {normalize_channel_name(channel): channel for channel in channels}

