# general imports
import numpy as np
from abc import ABC, abstractmethod

# local module imports
from . import layers as layers, constants
from . import global_functions as gf


class Membrane(ABC, layers.SolidLayer):
    def __new__(cls, membrane_dict, dx, **kwargs):
        model_type = membrane_dict.get('type', 'Constant')
        if model_type == 'Constant':
            return super(Membrane, cls).__new__(Constant)
        elif model_type == 'Linear':
            return super(Membrane, cls).__new__(LinearMembrane)
        elif model_type == 'Springer':
            return super(Membrane, cls).__new__(SpringerMembrane)
        elif model_type == 'YeWang2007':
            return super(Membrane, cls).__new__(YeWang2007Membrane)
        else:
            raise NotImplementedError('Specified membrane model not '
                                      'implemented. Available models are '
                                      'Constant, Linear, Springer, '
                                      'and YeWang2007.')

    def __init__(self, membrane_dict, dx, **kwargs):
        self.name = 'Membrane'
        membrane_dict['name'] = self.name
        super().__init__(membrane_dict, dx)

        # membrane temperature
        self.temp = np.zeros(self.dx.shape)

        # constant ionic conductivity of membrane
        self.ionic_conductivity = \
            membrane_dict.get('ionic_conductivity', 1.0e-3)

        # area specific membrane resistance
        self.omega_ca = np.zeros(self.dx.shape)

        # membrane resistance
        self.omega = np.zeros(self.dx.shape)

        self.calc_loss = membrane_dict.get('calc_loss', True)

        # voltage loss at the membrane
        self.v_loss = np.zeros(self.dx.shape)

        self.ionic_conductance = self.calc_conductance(self.ionic_conductivity)

        self.add_print_data(self.omega_ca, 'Membrane Resistance', 'Ohm-m²')

    @abstractmethod
    def calc_ionic_resistance(self, *args):
        pass

    def calc_voltage_loss(self, current_density, **kwargs):
        """
        Calculates the voltage loss at the membrane.
        """
        if not self.calc_loss:
            self.v_loss[:] = 0.
        else:
            self.v_loss[:] = self.omega_ca * current_density

    def update(self, current_density, humidity, *args):
        self.calc_ionic_resistance(humidity, *args)
        self.calc_voltage_loss(current_density)


class Constant(Membrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        # self.water_flux = np.zeros_like(self.dx)
        # water cross flux through the membrane
        self.omega[:] = 1.0 / self.ionic_conductance[0]
        self.omega_ca[:] = self.omega * self.area_dx

    def calc_ionic_resistance(self, *args):
        pass


class LinearMembrane(Membrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        self.basic_resistance = membrane_dict['basic_resistance']
        # basic electrical resistance of the membrane
        self.temp_coeff = membrane_dict['temperature_coefficient']
        # thermal related electrical resistance gain of the membrane

    def calc_ionic_resistance(self, *args):
        self.omega_ca[:] = \
            (self.basic_resistance - self.temp_coeff * self.temp)  # * 1e-2
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca


class WaterTransportMembrane(Membrane, ABC):

    FARADAY = constants.FARADAY

    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)

        # self.vapour_coeff = membrane_dict['vapour_transport_coefficient']
        # self.acid_group_conc = membrane_dict['acid_group_concentration']

        # Water cross flux through the membrane
        self.water_content = np.zeros((2, len(self.dx)))
        self.water_flux = np.zeros(self.dx.shape)
        self.diff_coeff = np.zeros(self.dx.shape)
        # Electro-osmotic drag coefficient
        self.eod = np.zeros(self.dx.shape)
        self.eod[:] = membrane_dict.get('electro-osmotic_drag_coeff', 1.07)

        self.add_print_data(self.water_flux,
                            'Membrane Water Flux', 'mol/(s-m²')

    @abstractmethod
    def calc_water_content(self, humidity):
        pass

    @abstractmethod
    def calc_diffusion_coefficient(self, *args):
        pass

    @abstractmethod
    def calc_cross_water_flux(self, current_density, humidity, *args):
        """
        Calculates the water cross flux through the membrane
        """
        pass

    def update(self, current_density, humidity, *args):
        self.calc_cross_water_flux(current_density, humidity)
        super().update(current_density, humidity, *args)


class SpringerMembrane(WaterTransportMembrane):
    """
    Models the ionic conductivity and water transport in the membrane
    (Nafion 117) according to the correlations proposed by:

    T.E. Springer, T.A. Zawodzinski, S. Gottesfeld. „Polymer Electrolyte Fuel
    Cell Model“. Journal of The Electrochemical Society 138, Nr. 8
    (1. August 1991): 2334–42. https://doi.org/10.1149/1.2085971.

    Currently, the water vapour activity at the membrane-catalyst interface
    is approximated by the channel humidity/water activity due to insufficient
    resolution of the through-plane concentration gradients.
    """
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        # Under-relaxation factor for water flux update
        self.urf = membrane_dict.get('underrelaxation_factor', 0.95)

    def calc_water_content(self, humidity):
        self.water_content[:] = \
            np.where(humidity < 1.0,
                     0.043 + 17.81 * humidity
                     - 39.85 * humidity ** 2. + 36. * humidity ** 3.0,
                     14.0 + 1.4 * (humidity - 1.0))
        return self.water_content

    def calc_diffusion_coefficient(self, *args):
        # Based on Springer et al. (1991);
        # as formulated by Kamarajugadda et al. (2008), non-continuous
        # function results in instabilities. Thus, the other formulation as
        # provided by Nguyen and White (1993) is used below.

        # Average water content for the diffusion coefficient calculation
        wc_avg = np.average(self.water_content, axis=0)
        # Minimum water content of the cathode and anode side
        wc_min = np.min(self.water_content, axis=0)
        diff_coeff_star = \
            np.where(wc_avg <= 2.0, 1.0,
                     np.where(wc_avg <= 3.0, 1.0 + 2.0 * (wc_avg - 2.0),
                              np.where(wc_avg <= 4.0,
                                       3.0 - 1.38 * (wc_avg - 3.0),
                                       2.563 - 0.33 * wc_avg
                                       + 0.0264 * wc_avg ** 2.0
                                       - 0.000671 * wc_avg ** 3.0)))

        diff_coeff = 1.0e-10 * np.exp(2416.0 * (1.0/303.0 * 1.0/self.temp)) \
            * diff_coeff_star

        # # Based on Nguyen and White (1993);
        # # as formulated by Kamarajugadda et al. (2008)
        # diff_coeff = 2.5/22.0 * 5.5e-11 * wc_avg \
        #     * np.exp(2416.0 * (1.0/303.0 * 1.0/self.temp))

        # Constant diffusion coefficient test
        # diff_coeff = 2.5e-6 * 0.0001

        self.diff_coeff[:] = diff_coeff
        return self.diff_coeff

    def calc_eod(self):
        # Minimum water content of the cathode and anode side
        wc_min = np.min(self.water_content, axis=0)
        # Electro-osmotic drag coefficient according to Springer et al. (1991),
        self.eod[:] = 2.5 * wc_min / 22.0
        return self.eod

    def calc_ionic_resistance(self, *args):
        """
        Calculates the membrane resistivity for Nafion 117
        according to Springer et al. (1991).
        """
        avg_water_content = np.average(self.water_content, axis=0)
        # water_content[water_content < 1.0] = 1.0
        # Membrane conductivity [S/m]
        mem_cond = (0.005139 * avg_water_content - 0.00326) \
            * np.exp(1268.0 * (0.0033 - 1. / self.temp)) * 1e2
        # Area-specific membrane resistance [Ohm-m²]
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        # Absolute resistance [Ohm]
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca

    def calc_cross_water_flux(self, current_density, humidity, *args):
        """
        Calculates the water cross flux through the membrane based on
        Springer et al. (1991). Water content on each membrane side is
        predicted by the channel humidities on each side, which probably
        overpredicts the difference significantly resulting in unreasonably
        high cross water fluxes and de-stabilizing the solution. For now, the
        humidities itself are used as the diffusion driving gradient,
        which should be reviewed further. More detailed resolution for
        calculation of humidities right at the membrane interface should
        resolve this issue in the future, which is part of a current project.
        """

        self.calc_water_content(humidity)
        self.calc_eod()
        self.calc_diffusion_coefficient()

        # Water flux due to electro-osmosis (ionic drag); -1 since
        # current_density is a magnitude but its direction should be
        # negative drag flux means from cathode to anode
        water_flux_drag = self.eod * -1.0 * current_density / self.FARADAY

        # Dry density of membrane (kg/m³), could be generalized as input
        rho_m = 2000.0
        # Molecular weight of membrane (kg/mol), could be generalized as input
        mw_m = 1.100

        # # Water flux due to diffusion as described by Springer et al. (1991),
        # # approximated as a linear gradient over the membrane thickness:
        # # water_content[0]: lambda at cathode side,
        # # water_content[1]: lambda at anode side
        # # negative diffusion flux mean from anode to cathode
        water_flux_diff = - rho_m / mw_m * self.diff_coeff \
            * (self.water_content[1] - self.water_content[0]) / self.thickness

        # Predicting the water content on each membrane side seems to be
        # overpredicting the gradient significantly, thus the humidities are
        # used for now
        # water_flux_diff = - rho_m / mw_m * diff_coeff \
        #     * (humidity[1] - humidity[0]) / self.thickness

        # Total water flux (diffusion based on temperature difference could be
        # added in the future)
        # Underrelaxation of flux
        water_flux = np.copy(self.water_flux)
        water_flux_new = water_flux_diff + water_flux_drag
        # water_flux_new = water_flux_drag
        self.water_flux[:] = self.urf * water_flux \
            + (1.0 - self.urf) * water_flux_new


class YeWang2007Membrane(SpringerMembrane):
    """
    Models the ionic conductivity and water transport in the membrane
    (GORE-Select) according to the correlations proposed by:

    X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
    Membrane-Electrode Assemblies: I. Membranes“. Journal of The Electrochemical
    Society 154, Nr. 7 (21. Mai 2007): B676. https://doi.org/10.1149/1.2737379.

    Currently, the water vapour activity at the membrane-catalyst interface
    is approximated by the channel humidity/water activity due to insufficient
    resolution of the through-plane concentration gradients.
    """
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        # Under-relaxation factor for water flux update
        self.urf = membrane_dict.get('underrelaxation_factor', 0.8)

    def calc_ionic_resistance(self, *args):
        """
        Equation 13 in:
        X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
        Membrane-Electrode Assemblies: I. Membranes“. Journal of The
        Electrochemical Society 154, Nr. 7 (21. Mai 2007): B676.
        https://doi.org/10.1149/1.2737379.
        """
        avg_water_content = np.average(self.water_content, axis=0)
        # water_content[water_content < 1.0] = 1.0
        # Membrane conductivity [S/m]
        mem_cond = 0.12 * avg_water_content ** 2.80 * 1e2
        # Area-specific membrane resistance [Ohm-m²]
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        # Absolute resistance [Ohm]
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca

    def calc_eod(self):
        self.eod[:] = 1.07
        return self.eod

    def calc_diffusion_coefficient(self, *args):
        """
        Equation 14 with k=0.5 in:
        X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
        Membrane-Electrode Assemblies: I. Membranes“. Journal of The
        Electrochemical Society 154, Nr. 7 (21. Mai 2007): B676.
        https://doi.org/10.1149/1.2737379.
        """
        super().calc_diffusion_coefficient()
        self.diff_coeff[:] *= 0.5
        return self.diff_coeff
