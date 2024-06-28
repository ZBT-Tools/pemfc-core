from .fluid import fluid as fl
from . import transport_layer as tl
from .fluid import diffusion_coefficient as dc
from . import discretization as dsct
from . import linear_system as ls


class DiffusionTransport:
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 fluid: (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 discretization: dsct.Discretization):

        self.dict = input_dict
        self.discretization = discretization
        self.fluid = fluid.copy(discretization.shape, plot_axis=-2)
        self.diff_coeff = dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas)
        self.transport_type = 'diffusion'
        self.transport_layers = [
            tl.TransportLayer.create(
                input_dict,
                {self.transport_type: [self.diff_coeff.d_eff[i],
                               self.diff_coeff.d_eff[i] * 1000,
                               self.diff_coeff.d_eff[i] * 1000]},
                discretization) for i in range(len(self.fluid.species_names))]

        self.linear_systems = [
            ls.BasicLinearSystem.create(item, self.transport_type) for item in
            self.transport_layers]
        print('test')

    # def update(self):
    #     # TODO: Define boundary and update input and update/solve diffusion
    #     #  problem