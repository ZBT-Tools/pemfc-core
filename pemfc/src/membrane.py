# General imports
import numpy as np
import scipy as sp
from abc import ABC, abstractmethod

# Local module imports
from . import transport_layer as tl, constants
# from pemfc.src import global_functions as gf
from . import discretization as dsct
from . import global_state as gs


class Membrane(tl.TransportLayer2D, ABC):
    # def __new__(cls, membrane_dict: dict, discretization: dsct.Discretization2D,
    #             *args, **kwargs):
    #     # model_type = membrane_dict.get('type', 'Constant')
    #     # transport_props = {
    #     #     'thermal': membrane_dict.get('thermal_conductivity', 0.0),
    #     #     'electrical': membrane_dict.get('electrical_conductivity', 0.0)}
    #     # if model_type == 'Constant':
    #     #     # return super(Membrane, cls).__new__(Constant, membrane_dict,
    #     #     #                                     transport_props,
    #     #     #                                     discretization)
    #     #     return super(Membrane, cls).__new__(Constant)
    #     #
    #     # elif model_type == 'Linear':
    #     #     return super(Membrane, cls).__new__(LinearMembrane)
    #     # elif model_type == 'Springer':
    #     #     return super(Membrane, cls).__new__(SpringerMembrane)
    #     # elif model_type == 'YeWang2007':
    #     #     return super(Membrane, cls).__new__(YeWang2007Membrane)
    #     # else:
    #     #     raise NotImplementedError('Specified membrane model not '
    #     #                               'implemented. Available models are '
    #     #                               'Constant, Linear, Springer, '
    #     #                               'and YeWang2007.')

    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        self.name = 'Membrane'
        membrane_dict['name'] = self.name
        solid_transport_properties = {
            'thermal': membrane_dict.get('thermal_conductivity', 0.0),
            'electrical': membrane_dict.get('electrical_conductivity', 0.0)}
        super().__init__(membrane_dict, solid_transport_properties,
                         discretization)

        # Membrane temperature
        self.temp = np.zeros(self.discretization.shape)

        # Constant ionic conductivity of membrane
        self.ionic_conductivity = \
            membrane_dict.get('ionic_conductivity', 1.0e-3)

        # area specific membrane resistance
        self.omega_ca = np.zeros(self.discretization.shape)

        # membrane resistance
        self.omega = np.zeros(self.discretization.shape)

        self.calc_loss = membrane_dict.get('calc_loss', True)

        # voltage loss at the membrane
        self.voltage_loss = np.zeros(self.discretization.shape)

        self.ionic_conductance = self.calc_conductance(self.ionic_conductivity)

        self.add_print_data(self.omega_ca, 'Membrane Resistance', 'Ohm-m²')

    @classmethod
    def create(cls, membrane_dict: dict, discretization: dsct.Discretization2D,
               *args, **kwargs):
        model_type = membrane_dict.get('type', 'Constant')
        if model_type == 'Constant':
            return Constant(membrane_dict, discretization)
        elif model_type == 'Linear':
            return LinearMembrane(membrane_dict, discretization)
        elif model_type == 'Springer':
            return SpringerMembrane(membrane_dict, discretization)
        elif model_type == 'YeWang2007':
            return YeWang2007Membrane(membrane_dict, discretization)
        elif model_type == 'VapourLiquidInterface':
            return VapourLiquidInterfaceResistanceMembrane(membrane_dict, discretization)
        else:
            raise NotImplementedError(
                'Specified membrane model is not implemented. Available models '
                'are Constant, Linear, Springer, and YeWang2007.')

    @abstractmethod
    def calc_ionic_resistance(self, *args, **kwargs):
        pass

    def calc_voltage_loss(self, current_density, **kwargs):

        """
        Calculates the voltage loss at the membrane.
        """
        if not self.calc_loss:
            self.voltage_loss[:] = 0.
        else:
            self.voltage_loss[:] = self.omega_ca * current_density

    def update(self, current_density, humidity, *args, **kwargs):
        self.calc_ionic_resistance(humidity, *args, **kwargs)
        self.calc_voltage_loss(current_density)


class Constant(Membrane):

    # def __new__(cls, membrane_dict: dict,
    #             discretization: dsct.Discretization2D, *args, **kwargs):
    #     instance = super().__new__(cls, membrane_dict, discretization,
    #                                *args, **kwargs)
    #     return instance

    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):

        super().__init__(membrane_dict, discretization, *args, **kwargs)
        # self.water_flux = np.zeros_like(self.dx)
        # water cross flux through the membrane
        self.omega[:] = 1.0 / self.ionic_conductance[0]
        self.omega_ca[:] = self.omega * self.discretization.d_area

    def calc_ionic_resistance(self, *args):
        pass


class LinearMembrane(Membrane):
    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        super().__init__(membrane_dict, discretization, *args, **kwargs)
        self.basic_resistance = membrane_dict['basic_resistance']
        # basic electrical resistance of the membrane
        self.temp_coeff = membrane_dict['temperature_coefficient']
        # thermal related electrical resistance gain of the membrane

    def calc_ionic_resistance(self, *args, **kwargs):
        self.omega_ca[:] = \
            (self.basic_resistance - self.temp_coeff * self.temp)  # * 1e-2
        self.omega[:] = self.omega_ca / self.discretization.d_area
        return self.omega, self.omega_ca


class WaterTransportMembrane(Membrane, ABC):

    FARADAY = constants.FARADAY

    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        super().__init__(membrane_dict, discretization, *args, **kwargs)

        # self.vapour_coeff = membrane_dict['vapour_transport_coefficient']
        # self.acid_group_conc = membrane_dict['acid_group_concentration']
        self.diff_coeff = np.zeros(self.discretization.shape)

        # Water cross flux through the membrane
        self.water_content = np.zeros((2, *self.discretization.shape))
        self.avg_water_content = np.zeros(self.discretization.shape)
        self.water_flux = np.zeros(self.discretization.shape)
        self.activity = np.zeros(self.water_content.shape)
        self.add_print_data(self.water_flux,
                            'Membrane Water Flux', 'mol/(s-m²)')
        self.add_print_data(self.avg_water_content,
                            'Membrane Water Content', '-')

    @abstractmethod
    def calc_water_content(self, humidity):
        pass

    @abstractmethod
    def calc_diffusion_coefficient(self, *args):
        pass

    @abstractmethod
    def calc_cross_water_flux(self, current_density, humidity, *args, **kwargs):
        """
        Calculates the water cross flux through the membrane
        """
        pass

    def update(self, current_density, humidity, *args, **kwargs):
        self.water_content[:] = self.calc_water_content(self.activity)
        self.avg_water_content[:] = np.average(self.water_content, axis=0)
        self.calc_cross_water_flux(current_density, humidity, *args, **kwargs)
        super().update(current_density, humidity, *args, **kwargs)


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
    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        super().__init__(membrane_dict, discretization, *args, **kwargs)
        # Under-relaxation factor for water flux update
        self.urf = membrane_dict.get('underrelaxation_factor', 0.95)

    def calc_water_content(self, activity):
        """
        Equations (16, 17) in
        Springer, T. E., T. A. Zawodzinski, S. Gottesfeld. „Polymer
        Electrolyte Fuel Cell Model“. Journal of The Electrochemical Society
        138, Nr. 8 (1. August 1991): 2334–42. https://doi.org/10.1149/1.2085971.
        """
        # water_content = \
        #     np.where(humidity < 1.0,
        #              0.043 + 17.81 * humidity
        #              - 39.85 * humidity ** 2. + 36. * humidity ** 3.0,
        #              14.0 + 1.4 * (humidity - 1.0))

        """
        Equation (25) in:
        Dickinson, E. J. F., G. Smith. „Modelling the Proton-Conductive 
        Membrane in Practical Polymer Electrolyte Membrane Fuel Cell (PEMFC) 
        Simulation: A Review“. Membranes 10, Nr. 11 (28. Oktober 2020): 310. 
        https://doi.org/10.3390/membranes10110310.
        """
        water_content = (
                0.3 + 6.0 * activity * (1.0 - np.tanh(activity - 0.5))
                + 3.9 * np.sqrt(activity)
                * (1.0 + np.tanh((activity - 0.89) / 0.23)))

        """
        linear test equation
        """
        # water_content = (humidity * 10.0)
        # Underrelaxation of water content
        return water_content

    def calc_diffusion_coefficient(self, *args):
        # Based on Springer et al. (1991);
        # as formulated by Kamarajugadda et al. (2008), non-continuous
        # function results in instabilities. Thus, the other formulation as
        # provided by Nguyen and White (1993) is used below.

        wc_avg = self.avg_water_content

        diff_coeff_star = (
            np.where(wc_avg <= 2.0, 1.0,
                     np.where(wc_avg <= 3.0, 1.0 + 2.0 * (wc_avg - 2.0),
                              np.where(wc_avg <= 4.0,
                                       3.0 - 1.38 * (wc_avg - 3.0),
                                       2.563 - 0.33 * wc_avg
                                       + 0.0264 * wc_avg ** 2.0
                                       - 0.000671 * wc_avg ** 3.0))))

        diff_coeff = (1.0e-10 * np.exp(2416.0 * (1.0 / 303.0 * 1.0 / self.temp))
                      * diff_coeff_star)

        # # Based on Nguyen and White (1993);
        # # as formulated by Kamarajugadda et al. (2008)
        # diff_coeff = 2.5/22.0 * 5.5e-11 * wc_avg \
        #     * np.exp(2416.0 * (1.0/303.0 * 1.0/self.temp))

        # Constant diffusion coefficient test
        # diff_coeff = 12e-10
        self.diff_coeff[:] = diff_coeff
        return diff_coeff

    def calc_eod(self):
        # Minimum water content of the cathode and anode side
        wc_min = np.min(self.water_content, axis=0)
        # Electro-osmotic drag coefficient according to Springer et al. (1991),
        eod = 2.5 * wc_min / 22.0
        return eod

    def calc_ionic_resistance(self, *args, **kwargs):
        """
        Calculates the membrane resistivity for Nafion 117
        according to Springer et al. (1991).
        """
        # water_content[water_content < 1.0] = 1.0
        # Membrane conductivity [S/m]

        mem_cond = (0.005139 * self.avg_water_content - 0.00326) \
            * np.exp(1268.0 * (0.0033 - 1. / self.temp)) * 1e2
        # Area-specific membrane resistance [Ohm-m²]
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        # Absolute resistance [Ohm]
        self.omega[:] = self.omega_ca / self.discretization.d_area
        return self.omega, self.omega_ca

    def calc_cross_water_flux(self, current_density, humidity, *args, **kwargs):
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

        # water_content = self.calc_water_content(humidity)
        # water_content_old = np.copy(self.water_content)
        # self.water_content[:] = (self.urf * water_content_old
        #                          + (1.0 - self.urf) * water_content)
        # self.avg_water_content[:] = np.average(self.water_content, axis=0)
        eod = self.calc_eod()
        diff_coeff = self.calc_diffusion_coefficient()

        # Water flux due to electro-osmosis (ionic drag); -1 since
        # current_density is a magnitude but its direction should be
        # negative corresponding to the direction from anode to cathode
        water_flux_drag = eod * -1.0 * current_density / self.FARADAY

        # Dry density of membrane (kg/m³), could be generalized as input
        rho_m = 2000.0
        # Molecular weight of membrane (kg/mol), could be generalized as input
        mw_m = 1.100

        # # Water flux due to diffusion as described by Springer et al. (1991),
        # # approximated as a linear gradient over the membrane thickness:
        # # water_content[0]: lambda at cathode side,
        # # water_content[1]: lambda at anode side
        # # negative diffusion flux mean from anode to cathode
        water_flux_diff = - rho_m / mw_m * diff_coeff \
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
        if gs.global_state.iteration == 150:
            pass
        self.water_flux[:] = self.urf * water_flux \
            + (1.0 - self.urf) * water_flux_new
        self.activity[:] = humidity
        return water_flux_new


class VapourLiquidInterfaceResistanceMembrane(SpringerMembrane):
    """
    Membrane model according to Springer, but with additional interface
    resistance for the water transport and water transport in the
    membrane according to the review in
    Dickinson, E. J. F., G. Smith. “Modelling the Proton-Conductive
    Membrane in Practical Polymer Electrolyte Membrane Fuel Cell (PEMFC)
    Simulation: A Review.” Membranes 10, no. 11 (October 28, 2020): 310.
    https://doi.org/10.3390/membranes10110310.
    """

    @staticmethod
    def calc_interface_saturation(liquid_pressure):
        """
        Equations (95) and (96) in:
        Dickinson, E. J. F., G. Smith. “Modelling the Proton-Conductive
        Membrane in Practical Polymer Electrolyte Membrane Fuel Cell (PEMFC)
        Simulation: A Review.” Membranes 10, no. 11 (October 28, 2020): 310.
        https://doi.org/10.3390/membranes10110310.
        """
        r_c = (4e-5 / liquid_pressure) * 1e9  # in nm
        saturation = 0.5 * (1.0 - sp.special.erf((np.log(r_c) - np.log(1.25))
                                                 / 0.3 * np.sqrt(2.0)))
        return saturation

    @staticmethod
    def calc_interfacial_resistance_coeff(humidity):
        """
        Equation (92)
        """
        # result = 4.4e-3
        return 1.04e-3 * np.exp(4.48e-4 * humidity)

    def calc_cross_water_flux(self, current_density, humidity, *args, **kwargs):
        """

        """
        liquid_pressure = kwargs.get('liquid_pressure', None)
        if liquid_pressure is None:
            liquid_saturation = 0.0
        else:
            liquid_saturation = self.calc_interface_saturation(liquid_pressure)

        liquid_equilibrium = kwargs.get('liquid_equilibrium', False)
        if liquid_equilibrium:
            liquid_saturation = 1.0

        # Use average liquid saturation due to lack of better knowledge
        # liquid_saturation_avg = np.average(liquid_saturation, axis=0)

        water_flux_new = super().calc_cross_water_flux(
            current_density, humidity)

        vapour_interfacial_resistance_coeff = (
            self.calc_interfacial_resistance_coeff(humidity))
        # Assumed factor of liquid-equilibrated interface resistance compared
        # to vapour-equilibrated interface resistance
        liquid_factor = 100.0
        vl_factor = (1.0 * (1.0 - liquid_saturation)
                     + liquid_saturation * liquid_factor)
        interfacial_resistance_coeff = (
                vl_factor * vapour_interfacial_resistance_coeff)


        # Recalculate membrane water activity
        activity = np.zeros(humidity.shape)
        activity[0] = (
            (interfacial_resistance_coeff[0] * humidity[0] - water_flux_new)
            / interfacial_resistance_coeff[0])
        activity[1] = (
            (interfacial_resistance_coeff[1] * humidity[1] + water_flux_new)
            / interfacial_resistance_coeff[1])
        activity = np.maximum(activity, 0.0)
        activity = np.minimum(activity, 1.0)

        #TODO: Missing link between activity and water content in equations
        self.activity[:] = activity
        # Underrelaxation of flux
        water_flux = np.copy(self.water_flux)
        # water_flux_new = water_flux_drag
        if gs.global_state.iteration == 150:
            pass
        self.water_flux[:] = self.urf * water_flux \
                             + (1.0 - self.urf) * water_flux_new

    def calc_eod(self):
        return 1.0

    def calc_ionic_resistance(self, *args, **kwargs):
        """
        Equation 13 in:
        X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
        Membrane-Electrode Assemblies: I. Membranes“. Journal of The
        Electrochemical Society 154, Nr. 7 (21. Mai 2007): B676.
        https://doi.org/10.1149/1.2737379.
        """
        avg_activity = np.average(self.activity, axis=0)
        # water_content[water_content < 1.0] = 1.0
        # Membrane conductivity [S/m]
        mem_cond = np.maximum(1e-1, 0.12 * avg_activity ** 2.80 * 1e2)
        # Area-specific membrane resistance [Ohm-m²]
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        # Absolute resistance [Ohm]
        self.omega[:] = self.omega_ca / self.discretization.d_area
        return self.omega, self.omega_ca

class YeWang2007Membrane(SpringerMembrane):
    """
    Models the ionic conductivity and water transport in the membrane
    (GORE-Select) according to the correlations proposed by:

    X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
    Membrane-Electrode Assemblies: I. Membranes“. Journal of The Electrochemical
    Society 154, Nr. 7 (21. Mai 2007): B676. https://doi.org/10.1149/1.2737379.
    """
    def __init__(self, membrane_dict: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        super().__init__(membrane_dict, discretization, *args, **kwargs)
        # Under-relaxation factor for water flux update
        self.urf = membrane_dict.get('underrelaxation_factor', 0.8)

    def calc_ionic_resistance(self, *args, **kwargs):
        """
        Equation 13 in:
        X. Ye, C.-Y. Wang. „Measurement of Water Transport Properties Through
        Membrane-Electrode Assemblies: I. Membranes“. Journal of The
        Electrochemical Society 154, Nr. 7 (21. Mai 2007): B676.
        https://doi.org/10.1149/1.2737379.
        """
        avg_activity = np.average(self.activity, axis=0)
        # water_content[water_content < 1.0] = 1.0
        # Membrane conductivity [S/m]
        mem_cond = np.maximum(1e-1, 0.12 * avg_activity ** 2.80 * 1e2)
        # Area-specific membrane resistance [Ohm-m²]
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        # Absolute resistance [Ohm]
        self.omega[:] = self.omega_ca / self.discretization.d_area
        return self.omega, self.omega_ca

    def calc_eod(self):
        return 1.0
        # return super().calc_eod()

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


# membrane_dict = {
#     'width': 3.0,
#     'length': 3.0,
#     'thickness': 0.1,
#     'electrical_conductivity': 1.0,
#     'thermal_conductivity': 1.0,
#     'ionic_conductivity': 15.0,
#     'type': 'Constant',
# }
# test_membrane = Membrane(membrane_dict, (20, 1))
# print(test_membrane.layer_dict)