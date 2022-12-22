# general imports
import numpy as np
from abc import ABC, abstractmethod

# local module imports
from . import layers as layers, constants


class Membrane(ABC, layers.SolidLayer):
    def __new__(cls, membrane_dict, dx, **kwargs):
        model_type = membrane_dict.get('type', 'Constant')
        if model_type == 'Constant':
            return super(Membrane, cls).__new__(Constant)
        elif model_type == 'Springer':
            return super(Membrane, cls).__new__(SpringerMembrane)
        elif model_type == 'Linear':
            return super(Membrane, cls).__new__(LinearMembrane)
        else:
            raise NotImplementedError('Specified membrane model not '
                                      'implemented. Available models are '
                                      'Constant, Springer, and Linear.')

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


class WaterTransportMembrane(Membrane, ABC):

    FARADAY = constants.FARADAY

    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)

        # self.vapour_coeff = membrane_dict['vapour_transport_coefficient']
        # self.acid_group_conc = membrane_dict['acid_group_concentration']
        # under-relaxation factor for water flux update
        self.urf = membrane_dict.get('underrelaxation_factor', 0.5)
        # water cross flux through the membrane
        self.water_content = np.zeros((2, len(self.dx)))
        self.water_flux = np.zeros(self.dx.shape)
        self.add_print_data(self.water_flux,
                            'Membrane Water Flux', 'mol/(s-m²')

    def calc_cross_water_flux(self, current_density, humidity):
        """
        Calculates the water cross flux through the membrane
        """
        pass

    def update(self, current_density, humidity, *args):
        self.calc_cross_water_flux(current_density, humidity)
        super().update(current_density, humidity, *args)


class Constant(Membrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        # self.water_flux = np.zeros_like(self.dx)
        # water cross flux through the membrane
        self.omega[:] = 1.0 / self.ionic_conductance[0]
        self.omega_ca[:] = self.omega * self.area_dx

    def calc_ionic_resistance(self, *args):
        pass


class SpringerMembrane(WaterTransportMembrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)

    def calc_ionic_resistance(self, humidity):
        """
        Calculates the membrane resistivity for Nafion 117
        according to (Springer, 1991).
        """
        humidity = np.average(humidity, axis=0)
        water_content = \
            np.where(humidity < 1.0,
                     0.043 + 17.81 * humidity
                     - 39.85 * humidity ** 2. + 36. * humidity ** 3.,
                     14.0 + 1.4 * (humidity - 1.0))
        water_content[water_content < 1.0] = 1.0
        mem_cond = (0.005139 * water_content - 0.00326) \
            * np.exp(1268.0 * (0.0033 - 1. / self.temp)) * 1e2
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca

    def calc_cross_water_flux(self, current_density, humidity):
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
        # humidity[humidity > 1.2] = 1.2
        water_content = \
            np.where(humidity < 1.0,
                     0.043 + 17.81 * humidity
                     - 39.85 * humidity ** 2. + 36. * humidity ** 3.0,
                     14.0 + 1.4 * (humidity - 1.0))

        self.water_content[:] = self.urf * self.water_content \
            + (1.0 - self.urf) * water_content


        # Average water content for the diffusion coefficient calculation
        wc_avg = np.average(water_content, axis=0)
        # Minimum water content of the cathode and anode side
        wc_min = np.min(water_content, axis=0)

        # Electro-osmotic drag coefficient according to Springer et al. (1991),
        eo_drag_coeff = 2.5 * wc_min / 22.0
        # # More recent publications found it to be close to unity (Xu et al.,
        # # 2017) or even decreasing with water content (Peng et al., 2011),
        # # therefore we simply use 1.0 for now (should be an input parameter)
        # eo_drag_coeff = 1.0

        # Water flux due to electro-osmosis (ionic drag); -1 since
        # current_density is a magnitude but its direction should be
        # negative drag flux means from cathode to anode
        water_flux_drag = eo_drag_coeff * -1.0 * current_density / self.FARADAY

        # Based on Springer et al. (1991);
        # as formulated by Kamarajugadda et al. (2008), non-continuous
        # function results in instabilities. Thus the other formulation as
        # provided by Nguyen and White (1993) is used below.
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
            * (water_content[1] - water_content[0]) / self.thickness

        # Predicting the water content on each membrane side seems to be
        # overpredicting the gradient significantly, thus the humidities are
        # used for now
        # water_flux_diff = - rho_m / mw_m * diff_coeff \
        #     * (humidity[1] - humidity[0]) / self.thickness

        # Total water flux (diffusion based on temperature difference could be
        # added in the future)
        # underrelaxation of flux
        water_flux = np.copy(self.water_flux)
        self.water_flux[:] = self.urf * water_flux \
            + (1.0 - self.urf) * (water_flux_diff + water_flux_drag)

        # # Previous calculation as presented by Chang et al. (2007)
        # dw = 2.1e-7 * np.exp(-2436. / self.temp)
        # water_content = 0.043 + 17.81 * humidity \
        #     - 39.85 * humidity ** 2. + 36. * humidity ** 3.
        #
        # divisor = 2. * self.vapour_coeff * self.acid_group_conc \
        #     * self.faraday_const
        # zeta_plus = \
        #     water_content[0] + water_content[1] + current_density / divisor
        # zeta_negative = (water_content[0] - water_content[1]
        #                  + 5. * current_density / divisor) \
        #     / (1. + dw * zeta_plus / (self.thickness * self.vapour_coeff))
        # m_c = 0.5 * (zeta_plus + zeta_negative)
        # m_a = 0.5 * (zeta_plus - zeta_negative)
        # self.w_cross_flow[:] = \
        #     current_density / self.faraday_const + self.acid_group_conc \
        #     * dw * (m_a ** 2. - m_c ** 2.) / (2. * self.thickness)


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
