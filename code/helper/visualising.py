import numpy as np
from bokeh.plotting import figure, show
from bokeh.layouts import column
from bokeh.palettes import Category10
from code.helper.utils import return_marker_channels


def plot_populations_flowkit(sample, populations, marker1, marker2, source="xform", marker_channels=None):
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
    channel1 = return_marker_channels(marker1, marker_channels)
    channel2 = return_marker_channels(marker2, marker_channels)

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