from abc import ABC, abstractmethod
import numpy as np

from .fluid import fluid as fl
from . import transport_layer as tl
from .fluid import diffusion_coefficient as dc
from . import discretization as dsct
from . import linear_system as ls
from . import global_functions as gf
from . import matrix_functions as mf
from . import half_cell as hc


class DiffusionTransport(ABC):
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 fluid: (fl.GasMixture, fl.CanteraGasMixture,
                         fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 discretization: dsct.Discretization, inert_id: int = -1):

        self.dict = input_dict
        self.discretization = discretization
        self.fluid = fluid.copy(discretization.shape, plot_axis=-2)
        if isinstance(fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            self.gas_diff_coeff = (
                dc.MixtureAveragedDiffusionCoefficient(self.fluid))
        elif isinstance(fluid, (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture)):
            self.gas_diff_coeff = (
                dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas))
        else:
            raise TypeError(
                'fluid must be of types (GasMixture, CanteraGasMixture, '
                'TwoPhaseMixture, or CanteraTwoPhaseMixture)')
        self.id_inert = inert_id
        self.transport_type = 'diffusion'
        input_dict['effective'] = True
        self.transport_layers = [
            tl.TransportLayer.create(
                input_dict,
                {self.transport_type: [self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i]]},
                discretization)
            for i in range(self.fluid.n_species) if i != self.id_inert]

        self.linear_systems = [
            ls.BasicLinearSystem.create(item, self.transport_type) for item in
            self.transport_layers]
        # test_lin_sys.set_neumann_boundary_conditions(15.0, axis=(0, 2),
        #                                              indices=(0, 1))
        # self.set_neumann_boundary_conditions(200.0, axis=(0, 2),
        #                                      indices=(0, 1))
        # self.set_dirichlet_boundary_conditions(15.0, axis=(0, 2),
        #                                        indices=(0, 0))

    @abstractmethod
    def update(self, *args, **kwargs):
        pass

    @classmethod
    def create(cls, input_dict: dict,
               fluid: (fl.GasMixture, fl.CanteraGasMixture,
                       fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
               discretization: dsct.Discretization, id_inert: int = -1):
        return GDLDiffusionTransport(input_dict, fluid, discretization,
                                     id_inert)


class GDLDiffusionTransport(DiffusionTransport):
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """

    def __init__(self, input_dict: dict,
                 fluid: (fl.GasMixture, fl.CanteraGasMixture,
                         fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 discretization: dsct.Discretization, inert_id: int = -1):
        super().__init__(input_dict, fluid, discretization, inert_id)
        self.initialize = True

        # Initialize Dirichlet boundary conditions at GDL-Channel interface
        # (z-index: 1); all other boundaries are flux boundaries;
        # specific boundary (rhs) values for GDL-Channel (concentration) and
        # CL-GDL interfaces (flux) will be set dynamically in update function
        for lin_sys in self.linear_systems:
            lin_sys.set_dirichlet_boundary_conditions(
                0.0, axes=(0, 2), indices=(0, 1))

    def update(self, temperature: np.ndarray, pressure: np.ndarray,
               molar_composition: np.ndarray, flux: np.ndarray):
        temperature = gf.rescale(temperature, self.fluid.temperature.shape)
        pressure = self.reshape_1d_input(pressure)
        pressure = gf.rescale(pressure, self.fluid.pressure.shape)

        molar_composition = np.asarray(
            [self.reshape_1d_input(mol_comp) for mol_comp in molar_composition])
        molar_composition = gf.rescale(
            molar_composition, self.fluid.gas.mole_fraction.shape)
        if self.initialize:
            self.fluid.update(temperature, pressure, molar_composition)
        else:
            self.fluid.update(temperature, pressure)

        # Update species transport for active species
        boundary_composition = [
            mol_comp for i, mol_comp in enumerate(molar_composition)
            if i != self.id_inert]
        boundary_flux = [flux[i] for i in range(len(flux)) if i !=
                         self.id_inert]
        self.gas_diff_coeff.update()
        for i, trans_layer in enumerate(self.transport_layers):
            trans_layer.update(
                {self.transport_type: [self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i]]})

        for i, lin_sys in enumerate(self.linear_systems):
            axes = (0,)
            indices = (1,)
            lin_sys.set_neumann_boundary_conditions(
                boundary_flux[i], axes=axes, indices=indices)
            axes = (0, 2)
            indices = (0, 1)
            concentration = mf.get_axis_values(
                    boundary_composition[i], axes=axes, indices=indices)
            lin_sys.set_dirichlet_boundary_conditions(
                concentration, axes=axes, indices=indices)
            lin_sys.update()

        # Update species transport for inert species
        if isinstance(self.fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            gas_constant = self.fluid.gas_constant
        else:
            gas_constant = self.fluid.gas.gas_constant
        total_gas_concentration = (
                self.fluid.pressure / (gas_constant * self.fluid.temperature))

        concentrations = [lin_sys.solution_array for lin_sys in self.linear_systems]
        active_concentrations = np.sum(concentrations, axis=0)
        inert_concentration = total_gas_concentration - active_concentrations
        concentrations.insert(self.id_inert, inert_concentration)
        concentrations = np.asarray(concentrations)
        self.initialize = False

    @staticmethod
    def reshape_1d_input(array):
        return np.moveaxis(np.asarray([[array,],]), (0, 1, 2), (0, 2, 1))
