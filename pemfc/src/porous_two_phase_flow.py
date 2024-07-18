import numpy as np
import json
import os
from porous_two_phase_flow import porous_layer as pl

from .fluid import (fluid as fl, diffusion_coefficient as dc,
                    evaporation_model as evap)
from . import transport_layer as tl
from . import discretization as dsct
from . import linear_system as ls
from . import global_functions as gf
from . import matrix_functions as mf
from . import global_state as gs
from . import diffusion_transport as dt


class TwoPhaseMixtureDiffusionTransport:
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 discretization: dsct.Discretization,
                 fluid: (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 *args, **kwargs):
        if not isinstance(fluid, fl.TwoPhaseMixture):
            raise TypeError('fluid must be of type TwoPhaseMixture')

        with open(os.path.join(gs.global_state.base_directory, 'settings',
                               'two_phase_settings.json')) as file:
            input_dict.update(json.load(file))
        input_dict['volume_fraction'] = (
            input_dict['porosity_model']['porosity'])
        input_dict['effective'] = True

        self.liquid_transport = dt.DiffusionTransport.create(
            input_dict, discretization, fluid)

        # self.fluid = fluid.copy(discretization.shape, plot_axis=-2)

        self.dict = input_dict
        # Create porosity model
        self.porosity_model = pl.PorousLayer(self.dict['porosity_model'],
                                             fluid)
        # Create liquid pressure diffusion transport system
        transport_property = (
                fluid.liquid.density / fluid.liquid.viscosity *
                np.asarray(self.dict['porosity_model']['permeability']))
        transport_layers = [tl.TransportLayer.create(
                input_dict, {'diffusion': transport_property}, discretization)]

    def update(self, *args, **kwargs):
        pass