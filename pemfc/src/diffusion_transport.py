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
        self.id_inert = inert_id
        self.transport_type = 'diffusion'
        self.fluid = fluid.copy(discretization.shape, plot_axis=-2)
        if 'diffusion_coefficient' in input_dict:
            input_dict['effective'] = False
            self.constant_diffusion = True
            self.transport_layers = [
                tl.TransportLayer.create(
                    input_dict,
                    {self.transport_type: input_dict['diffusion_coefficient']},
                    discretization)
                for i in range(self.fluid.n_species) if i != self.id_inert]
        else:
            self.constant_diffusion = False
            if isinstance(fluid, (fl.GasMixture, fl.CanteraGasMixture)):
                self.gas_diff_coeff = (
                    dc.MixtureAveragedDiffusionCoefficient(self.fluid))
            elif isinstance(fluid, (fl.TwoPhaseMixture,
                                    fl.CanteraTwoPhaseMixture)):
                self.gas_diff_coeff = (
                    dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas))
            else:
                raise TypeError(
                    'fluid must be of types (GasMixture, CanteraGasMixture, '
                    'TwoPhaseMixture, or CanteraTwoPhaseMixture)')
            input_dict['volume_fraction'] = input_dict['porosity']
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
        solution_shape = (len(self.linear_systems) + 1,
                          *self.linear_systems[0].solution_array.shape)
        self.concentrations = np.zeros(solution_shape)
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
        if 'boundary_patches' in input_dict:
            return GDLDiffusionTransport(input_dict, fluid, discretization,
                                         id_inert)
        else:
            raise NotImplementedError(
                'only GDLDiffusionTransport model is currently implemented '
                'which requires the subdict entry "boundary_patches" in the '
                '"input_dict" parameter')


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

        self.boundary_patches = input_dict['boundary_patches']
        neumann_axes = self.boundary_patches['Neumann']['axes']
        neumann_indices = self.boundary_patches['Neumann']['indices']
        neumann_shape = np.asarray(
            [mf.get_axis_values(conc, neumann_axes, neumann_indices)
             for i, conc in enumerate(self.concentrations)
             if i != self.id_inert]).shape
        self.flux_scaling_factors = np.zeros(neumann_shape)
        # Initialize Dirichlet boundary conditions at GDL-Channel interface
        # (z-index: 1); all other boundaries are flux boundaries;
        # specific boundary (rhs) values for GDL-Channel (concentration) and
        # CL-GDL interfaces (flux) will be set dynamically in update function
        # for lin_sys in self.linear_systems:
        #     lin_sys.set_dirichlet_boundary_conditions(
        #         0.0, axes=(0, 2), indices=(0, 1))

    def update(self, temperature: np.ndarray, pressure: np.ndarray,
               molar_composition: np.ndarray, flux: np.ndarray):
        temperature = gf.rescale(temperature, self.fluid.temperature.shape)
        pressure = self.reshape_1d_input(pressure)
        pressure = gf.rescale(pressure, self.fluid.pressure.shape)

        molar_composition = np.asarray(
            [self.reshape_input(mol_comp) for mol_comp in molar_composition])
        molar_composition = gf.rescale(
            molar_composition, self.fluid.gas.mole_fraction.shape)

        flux = np.asarray([self.reshape_input(item) for item in flux])
        flux = gf.rescale(flux, self.fluid.gas.mole_fraction.shape)
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
        if not self.constant_diffusion:
            self.gas_diff_coeff.update()
            for i, trans_layer in enumerate(self.transport_layers):
                trans_layer.update(
                    {self.transport_type: [self.gas_diff_coeff.d_eff[i],
                                           self.gas_diff_coeff.d_eff[i],
                                           self.gas_diff_coeff.d_eff[i]]})
            # trans_layer.conductance['diffusion'][:, 0:5, :, 0:10] = 1e-16
            # print(test)

        for i, lin_sys in enumerate(self.linear_systems):
            #
            # axes = (0,)
            # indices = (-1,)
            axes = self.boundary_patches['Neumann']['axes']
            indices = self.boundary_patches['Neumann']['indices']
            # lin_sys.set_neumann_boundary_conditions(
            #     boundary_flux[i], axes=axes, indices=indices)
            # axes = (0, 2)
            # indices = (0, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
            axes = self.boundary_patches['Dirichlet']['axes']
            indices = self.boundary_patches['Dirichlet']['indices']
            boundary_concentration = mf.get_axis_values(
                    boundary_composition[i], axes=axes, indices=indices)
            # lin_sys.set_dirichlet_boundary_conditions(
            #     boundary_concentration, axes=axes, indices=indices)
            boundary_conditions = self.boundary_patches
            boundary_conditions['Neumann']['values'] = boundary_flux[i]
            boundary_conditions['Dirichlet']['values'] = boundary_concentration
            lin_sys.update(rhs_values=boundary_conditions)

        # Update species transport for inert species
        if isinstance(self.fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            gas_constant = self.fluid.gas_constant
        else:
            gas_constant = self.fluid.gas.gas_constant
        total_gas_concentration = (
                self.fluid.pressure / (gas_constant * self.fluid.temperature))

        concentrations = [lin_sys.solution_array for
                          lin_sys in self.linear_systems]
        conc_array = np.array(concentrations)
        show_concentrations = np.round(np.moveaxis(conc_array,
                                       (0, 1, 2, 3), (0, 2, 1, 3)), 2)

        self.flux_scaling_factors[:] = self.calc_flux_scaling_factors(
            concentrations, boundary_composition)
        concentrations = [np.where(conc < 0.0, 0.0, conc) for conc in
                          concentrations]
        active_concentrations = np.sum(concentrations, axis=0)
        total_gas_concentration = gf.rescale(total_gas_concentration,
                                             active_concentrations.shape)
        inert_concentration = total_gas_concentration - active_concentrations
        concentrations.insert(self.id_inert, inert_concentration)
        concentrations = np.asarray(concentrations)
        self.concentrations[:] = concentrations
        show_concentrations = np.round(np.moveaxis(concentrations,
                                       (0, 1, 2, 3), (0, 2, 1, 3)), 2)
        conductances = [lin_sys.conductance for lin_sys in self.linear_systems]
        # show_conductances = np.round(np.)
        self.initialize = False

    def calc_flux_scaling_factors(
            self, neumann_bc_concentrations: (list, np.ndarray),
            dirichlet_bc_concentrations: (list, np.ndarray)):
        neumann_bc_axes = self.boundary_patches['Neumann']['axes']
        neumann_bc_indices = self.boundary_patches['Neumann']['indices']
        flux_boundary_concentrations = np.asarray([
            mf.get_axis_values(conc, neumann_bc_axes, neumann_bc_indices)
            for conc in neumann_bc_concentrations])
        channel_concentration = np.asarray([
            mf.get_axis_values(conc, neumann_bc_axes, neumann_bc_indices)
            for conc in dirichlet_bc_concentrations])
        flux_scaling_factors = (
            (channel_concentration - flux_boundary_concentrations)
            / channel_concentration)
        return flux_scaling_factors

    @staticmethod
    def reshape_1d_input(array: np.ndarray):
        return np.moveaxis(np.asarray([[array,],]), (0, 1, 2), (0, 2, 1))

    @staticmethod
    def reshape_2d_input(array: np.ndarray):
        return np.moveaxis(np.asarray([array,]), (0, 1, 2), (0, 1, 2))

    @classmethod
    def reshape_input(cls, array: np.ndarray):
        array = np.asarray(array)
        if array.ndim == 1:
            return cls.reshape_1d_input(array)
        elif array.ndim == 2:
            return cls.reshape_2d_input(array)
        else:
            raise ValueError('only 1- or 2-dimensional array allowed')
