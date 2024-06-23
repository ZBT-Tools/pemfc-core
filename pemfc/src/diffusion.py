from .fluid import fluid as fl
from .transport_layer import TransportLayer
from .fluid import diffusion_coefficient as dc
from .linear_system import BasicLinearSystem


class DiffusionTransport:
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 fluid: (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 transport_layer: TransportLayer):
        self.dict = input_dict
        self.shape = (2, *transport_layer.dsct.shape)
        self.fluid = fluid.copy(self.shape)
        self.transport_layer = transport_layer

        # for conc in self.fluid.gas.
        # self.lin_sys = Basi
        # Object of diffusion coefficient class to calculate diffusion
        # coefficients
        self.diff_coeff = dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas)
        # TODO: effective calculation of diff coefficient
        print('Initialized DiffusionTransport')
        self.diff_coeff_eff = self.diff_coeff


