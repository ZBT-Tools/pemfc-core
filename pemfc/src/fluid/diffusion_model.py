"""
Diffusion models
"""
import numpy as np
from abc import ABC, abstractmethod
from . import species
from . import fluid as fl
from pemfc.data import material_properties as mat_prop
from pemfc.src import constants as constants


class DiffusionModel(ABC):
    def __init__(self, fluid, knudsen=False, **kwargs):
        if not isinstance(fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            raise TypeError('Argument fluid of type GasMixture or '
                            'CanteraGasMixture must be provided')
        self.n_species = fluid.n_species
        self.species_names = fluid.species_names
        self.fluid = fluid
        array_shape = fluid.array_shape
        try:
            self.mw = {name: mat_prop.molecular_weight[name] for name in
                       self.species_names}
        except KeyError:
            raise NotImplementedError('Provided substance in species_list '
                                      'is not available (yet)')
        self.vd = {name: mat_prop.diffusion_volume[name] for name in
                   self.species_names}
        self.d_ij = np.zeros((self.n_species, self.n_species, *array_shape))
        self.knudsen = knudsen
        if knudsen:
            if 'pore_radius' in kwargs:
                self.pore_radius = kwargs.get('pore_radius')
            else:
                self.pore_radius = 1e-6
                print('Knudsen diffusion was activated, but argument '
                      'pore_radius not specified, default value of 1e-6 m is '
                      'used')
        self.updated = False

    def calc_binary_diffusion_coeff(self, species_i, species_j,
                                    temperature=None, pressure=None):
        """
        :return: binary diffusion coefficient for species_a and species_b (m²/s)
        :param species_i: string with sum formula for first species
        :param species_j: string with sum formula for second species
        :param temperature: temperature (K) (float or numpy array)
        :param pressure: pressure (Pa) (float or numpy array)

        According to:
        'Bao, Cheng, Zeyi Jiang, and Xinxin Zhang. “Modeling Mass Transfer in
        Solid Oxide Fuel Cell Anode: I. Comparison between Fickian,
        Stefan-Maxwell and Dusty-Gas Models.” Journal of Power Sources 310
        (April 2016): 32–40. https://doi.org/10.1016/j.jpowsour.2016.01.099.'
        """
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        const = 3.1976e-4
        t_factor = temperature ** 1.75
        mw_factor = np.sqrt(1.0 / self.mw[species_i] + 1.0 / self.mw[species_j])
        vd_factor = pressure * (self.vd[species_i] ** (1.0 / 3.0)
                                + self.vd[species_j] ** (1.0 / 3.0)) ** 2.0
        return const * t_factor * mw_factor / vd_factor

    def update_binary_diffusion_coeffs(self, temperature=None, pressure=None):
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        ij_list = []
        for i in range(self.n_species):
            for j in range(self.n_species):
                if [j, i] in ij_list:
                    self.d_ij[i, j] = self.d_ij[j, i]
                else:
                    d_ij = self.calc_binary_diffusion_coeff(
                        self.species_names[i], self.species_names[j],
                        temperature, pressure)
                    self.d_ij[i, j] = d_ij
                ij_list.append([i, j])

    def knudsen_diffusion_coeff(self, species_i, temperature=None):
        """
        :param species_i: string with sum formula for species
        :param temperature: temperature (K) (float or numpy array)
        :return: knudsen diffusion coefficient (m²/s)

        According to:
        'Bao, Cheng, Zeyi Jiang, and Xinxin Zhang. “Modeling Mass Transfer in
        Solid Oxide Fuel Cell Anode: I. Comparison between Fickian,
        Stefan-Maxwell and Dusty-Gas Models.” Journal of Power Sources 310
        (April 2016): 32–40. https://doi.org/10.1016/j.jpowsour.2016.01.099.'
        """
        if temperature is None:
            temperature = self.fluid.temperature
        if self.knudsen:
            return 2.0 / 3.0 * self.pore_radius \
                * np.sqrt(8.0 * constants.GAS_CONSTANT * temperature
                          / (np.pi * self.mw[species_i]))
        else:
            return None

    @abstractmethod
    def update(self, temperature=None, pressure=None, *args, **kwargs):
        # Calculate binary diffusion coefficients
        self.update_binary_diffusion_coeffs(temperature, pressure)
        self.updated = True


class MixtureAveragedDiffusionModel(DiffusionModel):
    def __init__(self, fluid):
        super().__init__(fluid)
        self.d_eff = np.zeros((self.n_species, *fluid.array_shape))

    def calc_diffusion_coeff(self, species_i, temperature=None, pressure=None,
                             mole_fractions=None, flux_ratio='stoich'):
        """
        :param species_i: string with sum formula for species
        :param temperature: temperature (K) (float or numpy array)
        :param pressure: pressure (Pa) (float or numpy array)
                :param mole_fractions: list or array of mole fractions for all
        species; for multiple dimensions, first dimension must be species;
        order must be equivalent to order in self.species list
        :param flux_ratio: ratio of molar fluxes between binary components,
        can be 'stoich' (equals unity), 'graham'

        According to:
        'Bao, Cheng, Zeyi Jiang, and Xinxin Zhang. “Modeling Mass Transfer in
        Solid Oxide Fuel Cell Anode: I. Comparison between Fickian,
        Stefan-Maxwell and Dusty-Gas Models.” Journal of Power Sources 310
        (April 2016): 32–40. https://doi.org/10.1016/j.jpowsour.2016.01.099.'
        """
        assert len(mole_fractions) == self.n_species
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        if mole_fractions is None:
            mole_fractions = self.fluid.mole_fraction
        x = np.asarray(mole_fractions)
        id_i = self.species_names.index(species_i)

        # Calculate binary diffusion coefficients
        if self.updated:
            d_ij = self.d_ij[id_i, :]
        else:
            d_ij = np.asarray(
                [self.calc_binary_diffusion_coeff(species_i, species_j,
                                                  temperature, pressure)
                    for species_j in self.species_names])
            self.d_ij[:] = d_ij

        if flux_ratio == 'stoich':
            beta_ij = np.ones(self.n_species) * -1
        elif flux_ratio == 'graham':
            beta_ij = np.asarray([np.sqrt(self.mw[species_i]/self.mw[species_j])
                                  for species_j in self.species_names]) * -1
        else:
            raise NotImplementedError

        inv_d = np.sum(np.asarray(
            [(x[id_j] - x[id_i] * beta_ij[id_j]) / d_ij[id_j]
             for id_j in range(self.n_species) if id_j != id_i]), axis=0)

        if self.knudsen:
            d_k_i = self.knudsen_diffusion_coeff(species_i, temperature)
            inv_d += 1.0 / d_k_i

        return 1.0 / inv_d

    def update(self, temperature=None, pressure=None, mole_fractions=None,
               flux_ratio='stoich', *args, **kwargs):
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        if mole_fractions is None:
            mole_fractions = self.fluid.mole_fraction
        super().update(temperature, pressure, *args, **kwargs)
        update_names = kwargs.get('update_names', None)
        if not isinstance(update_names, (tuple, list)):
            update_names = self.species_names
        # Calculate mixture averaged diffusion coefficients
        for i, name in enumerate(self.species_names):
            if name in update_names:
                self.d_eff[i] = \
                    self.calc_diffusion_coeff(name, temperature, pressure,
                                              mole_fractions,
                                              flux_ratio=flux_ratio)
        self.updated = True
