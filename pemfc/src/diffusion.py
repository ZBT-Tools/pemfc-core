from .fluid import fluid as fl
from .transport_layer import TransportLayer
from .fluid import diffusion_coefficient as dc
from . import discretization as dsct
from .linear_system import BasicLinearSystem


class DiffusionTransport:
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 fluid: (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 discretization: dsct.Discretization2D):

        self.dict = input_dict
        self.shape = (2, *discretization.shape)
        self.fluid = fluid.copy(self.shape, plot_axis=-2)


        self.diff_coeff = dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas)
        self.transport_layer = TransportLayer(
            input_dict, {'diffusion': [self.diff_coeff.d_eff,
                                       self.diff_coeff.d_eff * 1000,
                                       self.diff_coeff.d_eff * 1000]},
            discretization)
        print('test')
