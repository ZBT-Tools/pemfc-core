import sys, os
# https://stackoverflow.com/questions/16981921/relative-imports-in-python-3#comment122130933_16985066
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from  pemfc.src.fluid import fluid, species
from pemfc.src import channel as chl, flow_circuit, interpolation as ip, half_cell, constants
#reimporting channel because of channel_heat_transfer_simulation.py. Didn't want to cause breakage
from pemfc.src import channel
from pemfc.main_app import main

__all__ = ["fluid","species", "main", "chl", "flow_circuit", "ip", "half_cell","constants", "channel"]