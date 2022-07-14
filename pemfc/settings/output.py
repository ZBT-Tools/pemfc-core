""" Output options """
import os
file_path = os.path.dirname(__file__)

# output csv data (data output often takes the longest time)
save_csv_data = True

# output plots (data output is often takes the longest time)
save_plot_data = True

# show voltage losses in the U-I-graph (not working correctly at the moment;
# please remind to fix, if required)

directory = os.path.join(file_path, '..', '..', 'output')
