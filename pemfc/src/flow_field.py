# general imports
import numpy as np
from abc import ABC, abstractmethod

# local modul imports
from . import channel


class FlowField:
    """
    Class to calculate planar dimensions of flow field at the interface to the
    GDL. Given a channel object of class Channel, the outer dimensions 'width'
    and 'length' of the active area, the number of channels, 'channel_number',
    and optionally the 'rib_width', determine either the rib width or override
    the channel length to accommodate the given 'rib_width'
    """
    def __init__(self, name: str, flowfield_dict: dict, chl: channel.Channel):
        self.name = name
        self.channel = chl
        # number of channels of each half cell
        self.n_channels = flowfield_dict['channel_number']
        length = flowfield_dict['length']
        width = flowfield_dict['width']
        area = length * width

        # Rib width dictates and resets channel length if provided
        if 'rib_width' in flowfield_dict:
            self.rib_width = flowfield_dict['rib_width']
            self.channel.length = \
                area / ((chl.width + self.rib_width) * self.n_channels)
            area_factor = area / (chl.base_area * self.n_channels)
        else:
            area_factor = area / (chl.base_area * self.n_channels)
            if area_factor < 1.0:
                raise ValueError('width and length of cell result in a cell '
                                 'surface area  smaller than the area '
                                 'covered by channels')
            self.rib_width = chl.width * (area_factor - 1.0)

        self.width_straight_channels = \
            (chl.width + self.rib_width) * self.n_channels
        self.length_straight_channels = area / self.width_straight_channels
        self.active_area = area_factor * chl.base_area
        self.external_surface_factor = (width + length) \
            / (self.width_straight_channels + self.length_straight_channels)
        # self.active_area = area_factor * self.channel.base_area
        # factor active area with ribs / active channel area
        # self.active_area_dx = area_factor * chl.base_area_dx
