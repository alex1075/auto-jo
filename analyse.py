import bokeh
from bokeh.plotting import figure, output_file, show, column
import flowkit as fk
from fuzzywuzzy import fuzz
import numpy as np
from bokeh.palettes import Category10
from bokeh.models import Legend
from math import log10
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from gating import *


# Load the FCS file
sample = fk.Sample("fcs_files/CMV-SS_1.fcs")

# make channel list
channels = sample.pnn_labels
with open("channel_list.txt", "w") as f:
    for channel in channels:
        f.write(channel + "\n")

# Normalize channel names
def normalize_channel_name(channel):
    return channel.replace(" ", "").replace("*-A", "").replace("-A", "")

normalized_channels = {normalize_channel_name(channel): channel for channel in channels}




# match markers to channels in the panel
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

def return_marker_channels(marker):
    return marker_channels[marker]

def get_marker_index(marker):
    for i, channel in enumerate(channels):
        if marker in channel:
            return i
    return None

xform = fk.transforms.LogicleTransform(param_t=1024, param_w=0.5, param_m=4.5, param_a=0)
sample.apply_transform(xform)        




# Ensure the correct channels are mapped
channel1 = return_marker_channels("CD31")
channel2 = return_marker_channels("CD45")

# Select populations using the adjusted function
populations = get_marker_populations(sample, "CD31", "CD45", source="xform")

# Print the populations with f-string labels
for label, population in populations.items():
    print(f"{label}: {population.shape} events")

def plot_populations_flowkit(sample, populations, marker1, marker2, source="xform"):
    """
    Plots each population in FlowKit's scatter plot and creates an overarching plot.

    Args:
        sample (fk.Sample): The FlowKit sample object.
        populations (dict): A dictionary of populations as NumPy arrays, with f-string labels.
        marker1 (str): The name of the first marker (x-axis).
        marker2 (str): The name of the second marker (y-axis).
        source (str): Source of the data ('orig', 'raw', 'comp', or 'xform'). Default is 'xform'.
    """
    # Map marker names to channel names
    channel1 = return_marker_channels(marker1)
    channel2 = return_marker_channels(marker2)

    # Debug: Ensure channels are correctly mapped
    print(f"Using channels: {channel1} for {marker1}, {channel2} for {marker2}")

    # Get the full dataset
    full_data = sample.get_events(source=source)

    # Overarching plot
    all_plot = figure(
        title="Overarching Plot with All Populations",
        x_axis_label=marker1,
        y_axis_label=marker2,
        tools="pan,box_zoom,reset,save",
        width=800,
        height=800,
    )

    # Individual plots
    individual_plots = []
    colours = Category10[10]

    # Convert full_data and populations into structured arrays
    dtype = [(f"col{i}", full_data.dtype) for i in range(full_data.shape[1])]
    full_data_struct = full_data.view(dtype).reshape(full_data.shape[0])

    for i, (label, population_data) in enumerate(populations.items()):
        # Convert population data into structured array for comparison
        population_struct = population_data.view(dtype).reshape(population_data.shape[0])

        # Create a row-wise mask
        mask = np.isin(full_data_struct, population_struct)

        # Debug: Verify mask size
        print(f"Population: {label}, Mask Size: {mask.sum()}")

        # Individual plot for the population
        individual_plot = figure(
            title=f"Population: {label}",
            x_axis_label=marker1,
            y_axis_label=marker2,
            tools="pan,box_zoom,reset,save",
            width=400,
            height=400,
        )

        # Handle non-positive values if using log scale
        x_data = np.where(full_data[mask][:, 0] > 0, full_data[mask][:, 0], 1e-5)
        y_data = np.where(full_data[mask][:, 1] > 0, full_data[mask][:, 1], 1e-5)

        individual_plot.scatter(
            x=x_data,
            y=y_data,
            size=5,
            color=colours[i % len(colours)],
            alpha=0.6,
        )
        individual_plots.append(individual_plot)

        # Add population to the overarching plot
        all_plot.scatter(
            x=x_data,
            y=y_data,
            size=5,
            color=colours[i % len(colours)],
            alpha=0.6,
            legend_label=label,
        )

    # Customise the overarching plot legend
    all_plot.legend.title = "Populations"
    all_plot.legend.location = "top_right"
    all_plot.legend.click_policy = "hide"

    # Show all plots
    show(column(all_plot, *individual_plots))

# def plot_populations_flowkit(sample, populations, marker1, marker2, source="xform"):
#     """
#     Efficiently plots populations using FlowKit's sample.plot_scatter method.

#     Args:
#         sample (fk.Sample): The FlowKit sample object.
#         populations (dict): A dictionary of populations as NumPy arrays, with f-string labels.
#         marker1 (str): The name of the first marker (x-axis).
#         marker2 (str): The name of the second marker (y-axis).
#         source (str): Source of the data ('orig', 'raw', 'comp', or 'xform'). Default is 'xform'.
#     """
#     # Map marker names to channel names
#     channel1 = return_marker_channels(marker1)
#     channel2 = return_marker_channels(marker2)

#     # Debug: Ensure channels are correctly mapped
#     print(f"Using channels: {channel1} for {marker1}, {channel2} for {marker2}")

#     # Get the full dataset
#     full_data = sample.get_events(source=source)

#     # Convert full_data and population_data to structured arrays for efficient matching
#     dtype = [(f"col{i}", full_data.dtype) for i in range(full_data.shape[1])]
#     full_data_struct = full_data.view(dtype).reshape(full_data.shape[0])
    
#     for label, population_data in populations.items():
#         population_struct = population_data.view(dtype).reshape(population_data.shape[0])

#         # Create a mask using structured array matching
#         mask = np.isin(full_data_struct, population_struct)

#         # Debug: Verify the mask
#         print(f"Population: {label}, Mask Size: {mask.sum()}")

#         # Plot the population using FlowKit's scatter plot
#         sample.plot_scatter(
#             channel1,
#             channel2,
#             source=source,
#             event_mask=mask,
#         )

# Get populations
populations = get_marker_populations(sample, "CD31", "CD45", source="xform")

# Plot using Bokeh
plot_populations_flowkit(sample, populations, "CD31", "CD45", source="xform")

# Plot using FlowKit
# plot_populations_flowkit(sample, populations, "CD31", "CD45", source="xform")

output_file("CMV-SS_1.html")