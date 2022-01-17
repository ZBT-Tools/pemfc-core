# general imports
import numpy as np
from abc import ABC, abstractmethod

# local modul imports
from . import channel


class FlowField:
    def __init__(self, name, flowfield_dict, channel):
        self.name = name
        self.channel = channel
        # number of channels of each half cell
        self.n_channels = flowfield_dict['channel_number']
        length = flowfield_dict['length']
        width = flowfield_dict['width']
        area = length * width

        # rib width dictates and resets channel length if provided
        if 'rib_width' in flowfield_dict:
            self.rib_width = flowfield_dict['rib_width']
            self.width_straight_channels = \
                (channel.width + self.rib_width) * self.n_channels
            self.channel.length = \
                area / (channel.width + self.rib_width) * self.n_channels
            area_factor = area / (channel.base_area * self.n_channels)
        else:
            area_factor = area / (channel.base_area * self.n_channels)
            if area_factor < 1.0:
                raise ValueError('width and length of cell result in a cell '
                                 'surface area  smaller than the area '
                                 'covered by channels')
            self.rib_width = channel.width * (area_factor - 1.0)

        self.width_straight_channels = \
            (channel.width + self.rib_width) * self.n_channels
        self.length_straight_channels = area / self.width_straight_channels
        self.active_area = area_factor * channel.base_area
        # self.active_area = area_factor * self.channel.base_area
        # factor active area with ribs / active channel area
        self.active_area_dx = area_factor * channel.base_area_dx
        # self.active_area_dx = area_factor * self.channel.base_area_dx

