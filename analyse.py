from bokeh.plotting import output_file
import flowkit as fk
from code.helper.gating import *
from code.helper.utils import *
from code.helper.flow import *
from code.config.markers import markers
from code.helper.visualising import *


sample = read_sample("fcs_files/CMV-SS_1.fcs")

marker_channels = prepare_marker_channels(sample, markers, output_file="channel_list.txt")
xform = fk.transforms.LogicleTransform(param_t=1024, param_w=0.5, param_m=4.5, param_a=0)
sample.apply_transform(xform)        
print(marker_channels)
# Select populations using the adjusted function
populations = get_marker_populations(sample, "CD31", "CD45", source="raw", marker_channels=marker_channels)
# print(populations)
# Print the populations with f-string labels
# for label, population in populations.items():
    # print(f"{label}: {population.shape} events")

# Get populations
# populations = get_marker_populations(sample, "CD31", "CD45", source="xform", marker_channels=marker_channels)

# Plot using Bokeh
plot_populations_flowkit(sample, populations, "CD31", "CD45", source="xform", marker_channels=marker_channels)

output_file("CMV-SS_1.html")