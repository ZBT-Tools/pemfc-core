"""
Diffusion models
"""
import numpy as np
from abc import ABC, abstractmethod
from . import species
from . import fluid as fl
from ...data import material_properties as mat_prop
from ...src import constants as constants


class EvaporationModel(ABC):
    def __new__(cls, fluid, model_dict, *arg, **kwargs):
        model_type = model_dict.get('type', 'HertzKnudsenSchrage')
        if model_type == 'HertzKnudsenSchrage':
            return super(EvaporationModel, cls).\
                __new__(HertzKnudsenSchrageModel)
        elif model_type == 'WangSi':
            return super(EvaporationModel, cls).__new__(WangSiModel)
        else:
            raise NotImplementedError

    def __init__(self, fluid, model_dict, *arg, **kwargs):
        if not isinstance(fluid, (fl.TwoPhaseMixture,
                                  fl.CanteraTwoPhaseMixture)):
            raise TypeError('Argument fluid of type TwoPhaseMixture or '
                            'CanteraTwoPhaseMixture must be provided')
        self.dict = model_dict
        self.name = model_dict.get('name', None)
        self.n_species = fluid.n_species
        self.species_names = fluid.species_names
        self.fluid = fluid
        self.mw = fluid.gas.species.mw[fluid.id_pc]
        self.updated = False

    @abstractmethod
    def calc_evaporation_rate(self, temperature=None, pressure=None,
                              *args, **kwargs):
        """
        :param temperature: temperature (K) (float or numpy array)
        :param pressure: pressure (Pa) (float or numpy array)
        """
        return None

    @abstractmethod
    def update(self, temperature=None, pressure=None, *args, **kwargs):
        # Calculate binary diffusion coefficients
        self.updated = True


class HertzKnudsenSchrageModel(EvaporationModel):
    def __init__(self, fluid, model_dict, *arg, **kwargs):
        super().__init__(fluid, model_dict)
        self.evap_coeff = model_dict['evaporation_coefficient']
        self.mw = fluid.gas.species.mw[fluid.id_pc]
        self.evap_rate = np.zeros(*fluid.array_shape)
        self. evap_coeff_factor = 2.0 * self.evap_coeff \
            / (2.0 - self.evap_coeff)

        k_b = constants.BOLTZMANN
        m = self.mw / constants.AVOGADRO_NUMBER
        self.const_factor = np.sqrt(m / (2.0 * np.pi * k_b))

    def calc_evaporation_rate(self, temperature=None, pressure=None,
                              capillary_pressure=None,
                              temperature_liquid=None, *args, **kwargs):
        """
        :param temperature: temperature (K) (float or numpy array)
        :param pressure: pressure (Pa) (float or numpy array)
        :param capillary_pressure: capillary_pressure (Pa) (float or numpy
        array), if given is used to correct the saturation pressure
        :param temperature_liquid: liquid phase temperature (K), if not provided
        equilibrium conditions are assumed and the temperature jump across
        the vapor and liquid interface is neglected
        :rtype: evaporation/condensation rate (kg / (m²-s)) (float or numpy
        array)

        According to:
        'Safi, Mohammad Amin, Nikolaos I. Prasianakis, John Mantzaras,
        Adrien Lamibrac, and Felix N. Büchi. “Experimental and Pore-Level
        Numerical Investigation of Water Evaporation in Gas Diffusion Layers of
        Polymer Electrolyte Fuel Cells.” International Journal of Heat and Mass
        Transfer 115 (December 2017): 238–49.
        https://doi.org/10.1016/j.ijheatmasstransfer.2017.07.050.'
        """
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure

        if temperature_liquid is not None:
            temperature_liq = temperature_liquid
            temperature_vap = temperature
        else:
            temperature_liq = temperature
            temperature_vap = temperature

        p_sat = self.fluid.saturation_pressure
        if capillary_pressure is not None:
            gas_const = constants.GAS_CONSTANT
            p_sat *= \
                np.exp(capillary_pressure * self.mw
                       / (gas_const * temperature * self.fluid.gas.density))
        p_vap = pressure * self.fluid.mole_fraction[self.fluid.id_pc]
        return self.evap_coeff_factor * self.const_factor \
            * (p_sat / np.sqrt(temperature_liq)
               - p_vap / np.sqrt(temperature_vap))

    def update(self, temperature=None, pressure=None, *args, **kwargs):
        super().update(temperature, pressure, *args, **kwargs)
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        self.evap_rate[:] = \
            self.calc_evaporation_rate(temperature, pressure, *args, **kwargs)


class WangSiModel(EvaporationModel):
    def __init__(self, fluid, model_dict, *arg, **kwargs):
        super().__init__(fluid, model_dict)
        self.evap_coeff = model_dict['evaporation_coefficient']
        self.cond_coeff = model_dict['condensation_coefficient']
        self.evap_rate = np.zeros(*fluid.array_shape)
        self. evap_coeff_factor = 2.0 * self.evap_coeff \
            / (2.0 - self.evap_coeff)

        # k_b = constants.BOLTZMANN
        # m = self.mw / constants.AVOGADRO_NUMBER
        self.gas_constant = constants.GAS_CONSTANT
        # self.const_factor = np.sqrt(m / (2.0 * np.pi * k_b))

    def calc_evaporation_rate(self, saturation=None, temperature=None,
                              pressure=None, capillary_pressure=None,
                              *args, **kwargs):
        """
        :param saturation: liquid water saturation (-) (float or numpy array)
        :param temperature: temperature (K) (float or numpy array)
        :param pressure: pressure (Pa) (float or numpy array)
        :param capillary_pressure: capillary_pressure (Pa) (float or numpy
        array), if given is used to correct the saturation pressure
        :rtype: evaporation/condensation rate (kg / s) (float or numpy array)

        According to:
        'Si, Chao, Gui Lu, Xiao-Dong Wang, and Duu-Jong Lee.
        “Gas Diffusion Layer Properties on the Performance of Proton Exchange
        Membrane Fuel Cell: Pc-s Relationship with K-Function.” International
        Journal of Hydrogen Energy 41, no. 46 (December 2016): 21827–37.
        https://doi.org/10.1016/j.ijhydene.2016.07.005.
        """
        if saturation is None:
            raise ValueError('saturation must be provided')
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure

        if 'porosity' in kwargs:
            porosity = kwargs['porosity']
        else:
            raise ValueError('porosity must be provided in kwargs')

        p_sat = self.fluid.saturation_pressure
        if capillary_pressure is not None:
            gas_const = constants.GAS_CONSTANT
            p_sat *= \
                np.exp(capillary_pressure * self.mw
                       / (gas_const * temperature * self.fluid.gas.density))
        p_vap = pressure * self.fluid.mole_fraction[self.fluid.id_pc]

        p_diff = p_vap - p_sat
        with np.errstate(invalid='ignore'):
            pressure_ratio = \
                np.where(p_diff == 0.0, 1.0, np.abs(p_diff) / p_diff)

        mw_h2o = self.fluid.gas.species.mw[self.fluid.id_pc]
        x_h2o = self.fluid.mole_fraction[self.fluid.id_pc]

        h_pc_c = self.cond_coeff * porosity * (1.0 - saturation) * \
            x_h2o / (2.0 * self.gas_constant * temperature) * \
            (1.0 + pressure_ratio)
        h_pc_e = self.evap_coeff * porosity * saturation * \
            self.fluid.liquid.density / (2.0 * mw_h2o) * \
            (1.0 - pressure_ratio)
        return -(h_pc_c + h_pc_e) * mw_h2o * p_diff

    def update(self, temperature=None, pressure=None, *args, **kwargs):
        super().update(temperature, pressure, *args, **kwargs)
        if temperature is None:
            temperature = self.fluid.temperature
        if pressure is None:
            pressure = self.fluid.pressure
        self.evap_rate[:] = \
            self.calc_evaporation_rate(temperature, pressure, *args, **kwargs)
