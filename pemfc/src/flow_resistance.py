# general imports
from abc import ABC, abstractmethod
import numpy as np

# local module imports
from . import interpolation as ip
from . import constants


class FlowResistance(ABC):
    """
    Abstract base class for flow resistances creating pressure drop in channels.
    A referenced flow channel object must be injected at instantiation.
    Concrete-class creation is based on the "type" parameter defined in the
    provided dictionary "zeta_dict".
    """
    def __new__(cls, channel, zeta_dict, **kwargs):
        zeta_type = zeta_dict.get('type', 'Constant')
        if zeta_type == 'Constant':
            return super(FlowResistance, cls).\
                __new__(ConstantFlowResistance)
        elif zeta_type == 'WallFriction':
            return super(FlowResistance, cls).\
                __new__(WallFrictionFlowResistance)
        elif zeta_type == 'RennelsTeeMain':
            return super(FlowResistance, cls).\
                __new__(RennelsTeeMainFlowResistance)
        elif zeta_type == 'RennelsTeeBranch':
            return super(FlowResistance, cls).\
                __new__(RennelsTeeBranchFlowResistance)
        else:
            raise NotImplementedError

    def __init__(self, channel, zeta_dict, **kwargs):
        self.channel = channel
        self.value = None

    def update(self):
        """
        Updates resistance values based on local flow conditions
        """
        pass

    @abstractmethod
    def calc_pressure_drop(self, **kwargs):
        """
        Calculate element-wise pressure drop depending on the resistance values
        :return:
        """
        pass


class ConstantFlowResistance(FlowResistance):
    """
    Simple class for a flow resistance based on a constant zeta
    value. The resistance will be considered for each discrete element. So
    single resistances for a channel can only be considered by averaging them
    over all elements with this class.
    """
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.value = zeta_dict['value']

    def calc_pressure_drop(self, reference_velocity='average'):
        """
        Calculates element-wise pressure drop with constant resistance value
        (zeta value). Parameter "reference_velocity" determines the used
        velocity. Only important if velocity variations occur along channels.
        :param reference_velocity: average (default), in, out
        :return: element-wise pressure drop due to constant friction resistance
        """
        # Node-based pressure drop
        dp_node = 0.5 * self.channel.fluid.density * self.value \
            * self.channel.velocity ** 2.0

        if reference_velocity == 'in':
            ref_factor = 1
        elif reference_velocity == 'out':
            ref_factor = -1
        else:
            ref_factor = 0

        slice_direction = int(ref_factor * self.channel.flow_direction)
        if slice_direction > 0:
            return dp_node[:-1]
        elif slice_direction < 0:
            return dp_node[1:]
        else:
            return ip.interpolate_1d(dp_node)


class WallFrictionFlowResistance(FlowResistance):
    """
    Class for a flow resistance based on wall friction. Darcy-Weisbach
    coefficient is calculated for each element based on the local flow
    conditions ranging from laminar, transitional to turbulent.
    value. The Blasius equation is used for calculation of the friction factor
    at the moment, but could be easily extended.
    """
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.method = zeta_dict.get('method', 'Blasius')
        self.value = np.zeros(self.channel.n_nodes)

    def update(self):
        """
        Updates resistance values based on local flow conditions
        :return: None
        """
        # reynolds = ip.interpolate_1d(self.channel.reynolds)
        reynolds = self.channel.reynolds
        lam = reynolds < 2200.0
        turb = np.invert(lam)
        lam_id = np.where(lam)
        turb_id = np.where(turb)

        # Darcy friction factor
        factor = np.zeros(reynolds.shape)
        if np.any(lam):
            reynolds_lam = reynolds[lam_id]
            if self.channel.aspect_ratio == 1.0:
                f_reynolds = 64.0
            elif self.channel.cross_shape == 'rectangular':
                eps = self.channel.aspect_ratio
                f_reynolds = 4.0 * 24.0 \
                    / ((1.0 + eps) ** 2.0
                       * (1.0 - (192.0 * eps / np.pi ** 5.0
                                 * np.tanh(np.pi / (2.0 * eps)))))
            elif self.channel.cross_shape in ('triangular', 'trapezoidal'):
                f_reynolds = 64.0
            else:
                raise NotImplementedError

            factor_lam = \
                np.divide(f_reynolds, reynolds_lam,
                          where=reynolds_lam > 1e-6)
            factor[lam_id] = factor_lam
        if np.any(turb):
            reynolds_turb = reynolds[turb_id]
            if self.method == 'Blasius':
                factor[turb_id] = 0.3164 * np.power(reynolds_turb, -0.25)
            else:
                raise NotImplementedError
        np.seterr(under='ignore')
        self.value[:] = self.channel.dx_node / self.channel.d_h * factor
        self.value[np.isnan(self.value)] = 0.0
        self.value[self.value < constants.SMALL] = 0.0
        np.seterr(under='raise')
        # if np.any(self.value > 1e3):
        #     raise FloatingPointError

    def calc_pressure_drop(self):
        """
        Calculates element-wise pressure drop based on Darcy friction factor.
        Nodal-based resistance coefficients are rescaled to element-based
        resistance coefficients with their respective diameters and lengths.
        Relevant for channels with varying diameters and fluid properties.
        :return: element-wise pressure drop due to wall friction
        """
        dp_node = 0.5 * self.channel.fluid.density * self.value \
                  * self.channel.velocity ** 2.0
        # Weighting according to dx and dx_node lengths for element-wise
        # pressure drop
        dp_node_1 = dp_node[:-1]
        dp_node_2 = dp_node[1:]
        dx_node_1 = self.channel.dx_node[:-1]
        dx_node_2 = self.channel.dx_node[1:]
        dx = self.channel.dx
        # dp = (dx_node_1 * dp_node_1 + dx_node_2 * dp_node_2) / dx
        dp = (dp_node_1 + dp_node_2) * dx / (dx_node_1 + dx_node_2)
        return dp


class JunctionFlowResistance(FlowResistance, ABC):
    """
    Abstract base class to calculate element-wise additional pressure
    drop due to flow branching off and combining at T-junctions located at each
    element. The resistance calculated here only applies for the pressure
    drop within the main channel, not from the channel to the branches.
    Reference location is always at the higher velocity in the main
    channel, i.e. the combined flow velocity.
    Concrete implementations of this class must define the method
    "calc_resistance_values"

    :return: element-wise additional pressure drop due to flow branching off or
    combining via T-junctions
    """
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        # self.coeffs = zeta_dict['coefficients']
        self.value = np.zeros(self.channel.n_ele)
        self.main_flow_ratio = np.zeros(self.channel.n_ele)
        self.velocity = np.zeros(self.channel.n_ele)
        self.density = np.zeros(self.channel.n_ele)
        self.mass_source_sign = np.zeros(self.channel.n_ele)
        self.pressure = np.zeros(self.channel.n_ele)

    def update(self):
        """
        Updates resistance values and determine reference velocity based on
        local flow conditions. Always the combined velocity is used to
        calculate the pressure drop. For dividing T-junctions this is the
        inflow velocity, for combining T-junctions this is the outflow
        velocity. For single channel object with this resistance model both
        types are handled simultaneously.
        :return: None
        """

        # Vectorized method to determine the reference velocities and branch
        # flow direction (dividing or combining) at each
        # T-junction at each element
        if hasattr(self.channel.fluid, 'species'):
            mass_source = np.sum(self.channel.mass_source, axis=0)
        else:
            mass_source = self.channel.mass_source

        nodes = self.channel.n_nodes
        flow_direction = self.channel.flow_direction
        self.mass_source_sign[:] = np.sign(mass_source)
        # List of indices in vector which determine the combined flow
        combined_ids = np.arange(1, nodes) \
            - np.where(np.int64(self.mass_source_sign * flow_direction) <= 0,
                       1, 0)
        # List of indices in vector which determine the reduced flow
        reduced_ids = np.arange(1, nodes) \
            - np.where(np.int64(self.mass_source_sign * -flow_direction) <= 0,
                       1, 0)

        # Resistance values are defined with respect to the combined flow point
        self.velocity[:] = self.channel.velocity[combined_ids]
        density = self.channel.fluid.density
        self.density[:] = density[combined_ids]
        self.pressure[:] = self.channel.pressure[combined_ids]
        mass_flow = self.channel.mass_flow_total
        combined_volume_flow = mass_flow[combined_ids] * density[combined_ids]
        reduced_volume_flow = mass_flow[reduced_ids] * density[reduced_ids]
        self.main_flow_ratio[:] = \
            np.divide(reduced_volume_flow, combined_volume_flow,
                      out=self.main_flow_ratio,
                      where=combined_volume_flow != 0.0)

        # np.polynomial.polynomial.polyval(self.main_flow_ratio,
        #                                  self.coeffs)

        self.value[:] = self.calc_resistance_values()

    @abstractmethod
    def calc_resistance_values(self, **kwargs):
        """
        Calculates resistance values based on given correlation. Must be
        defined in concrete class implementation.
        :param kwargs:
        :return: vector including the calculated element-wise resistance values
        along the channel
        """
        pass

    def calc_pressure_drop(self):
        """
        Calculates additional element-wise pressure drop due to T-junctions at
        each element with flow dividing into branches or mixing with flow
        coming from branches. Resistance values must be updated previously
        according to local flow conditions with update-method.

        :return: additional element-wise pressure drop due to T-Junctions
        """
        return 0.5 * self.density * self.value * self.velocity ** 2.0


class RennelsTeeMainFlowResistance(JunctionFlowResistance):
    """
    Concrete implementation of the abstract JunctionFlowResistance class  
    defining the "calc_resistance_values" with values from the publication:
    Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
    Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
    https://doi.org/10.1002/9781118275276.
    Equations (16.5) for dividing flow and Equation (16.22) for combining
    flow to calculate the resistance values in the main run are implemented.
    The branch diameter must be provided as the parameter "branch_diameter"
    in the zeta_dict dictionary. Optionally, a fitting radius with the
    parameter "fitting_radius" can be supplied. The default value is zero.
    """
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.branch_diameter = zeta_dict["branch_diameter"]
        self.fitting_radius = zeta_dict.get("fitting_radius", 0.0)
        r_to_d = self.fitting_radius / self.branch_diameter
        self.C_M = 0.23 + 1.46 * r_to_d - 2.75 * r_to_d ** 2.0 \
            + 1.65 * r_to_d ** 3.0
        self.C_xC = 0.08 + 0.56 * r_to_d - 1.75 * r_to_d ** 2.0 \
            + 1.83 * r_to_d ** 2.0

    @staticmethod
    def calc_dividing_junction_value(flow_ratio):
        """
        Equation (16.5) from
        Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
        Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
        https://doi.org/10.1002/9781118275276.
        :param flow_ratio: volumetric flow ratio array in main run
        (w2 / w1 in publication)
        :return: zeta resistance value
        """
        return 0.36 - 0.98 * flow_ratio + 0.62 * flow_ratio ** 2.0 \
            + 0.03 * flow_ratio ** 8

    def calc_combining_junction_value(self, flow_ratio):
        """
        Equation (16.22) from
        Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
        Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
        https://doi.org/10.1002/9781118275276.
        :param flow_ratio: volumetric flow ratio array in main run
        (w2 / w1 in publication)
        :return: zeta resistance value
        """
        return 1.0 - 0.95 * flow_ratio ** 2.0 \
            - 2.0 * self.C_xC * (flow_ratio - flow_ratio ** 2.0) \
            - 2.0 * self.C_M * (1.0 - flow_ratio)

    def calc_resistance_values(self, **kwargs):
        """
        Calculates resistance values based on given correlation. Must be
        defined in concrete class implementation.
        :param kwargs:
        :return: vector including the calculated element-wise resistance values
        along the channel
        """
        # Calculate both resistance values, simpler due to using the
        # parameters self.C_xC and self.C_M which could be arrays as well. If
        # only the values at the required positions should be calculated,
        # the parameters should be made function inputs as well.
        value = np.zeros(self.value.shape)
        value_dividing = \
            self.calc_dividing_junction_value(self.main_flow_ratio)
        value_combining = \
            self.calc_combining_junction_value(self.main_flow_ratio)
        # Determine indices for dividing or combining elements
        id_dividing = np.nonzero(self.mass_source_sign <= 0)
        id_combining = np.nonzero(self.mass_source_sign > 0)
        # Assign result values to corresponding elements
        value[id_dividing] = value_dividing[id_dividing]
        value[id_combining] = value_combining[id_combining]
        return value


class RennelsTeeBranchFlowResistance(RennelsTeeMainFlowResistance):
    """
    Concrete implementation of the abstract JunctionFlowResistance class
    defining the "calc_resistance_values" with values from the publication:
    Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
    Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
    https://doi.org/10.1002/9781118275276.
    Equations (16.13) for dividing flow and Equation (16.30) for combining
    flow to calculate the resistance values for the additional pressure drop
    from the main flow to the branch are implemented.
    The branch diameter must be provided as the parameter "branch_diameter"
    in the zeta_dict dictionary. Optionally, a fitting radius with the
    parameter "fitting_radius" can be supplied. The default value is zero.
    """
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.diameter_ratio = self.branch_diameter / self.channel.diameter
        r_to_d = self.fitting_radius / self.branch_diameter
        self.K_9_3 = 0.57 - 1.07 * r_to_d ** 0.5 - 2.13 * r_to_d \
            + 8.24 * r_to_d ** 1.5 - 8.48 * r_to_d ** 2.0 \
            + 2.9 * r_to_d ** 2.5
        self.C_yC = 1.0 - 0.25 * self.diameter_ratio ** 1.3 \
            - (0.11 * r_to_d - 0.65 * r_to_d ** 2.0
               + 0.83 * r_to_d ** 3.0) * self.diameter_ratio ** 2.0

    def calc_dividing_junction_value(self, flow_ratio):
        """
        Equation (16.13) from
        Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
        Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
        https://doi.org/10.1002/9781118275276.
        :param flow_ratio: volumetric flow ratio array for branch flow to
        combined main flow (w3 / w1 in publication)
        :return: zeta resistance value
        """
        return 1.0 - 1.13 * flow_ratio \
            + (0.81 + (1.12 * self.diameter_ratio
                       - 1.08 * self.diameter_ratio ** 3 + self.K_9_3) *
               self.diameter_ratio ** -4.0) * flow_ratio ** 2.0

    def calc_combining_junction_value(self, flow_ratio):
        """
        Equation (16.30) from
        Rennels, Donald C., and Hobart M. Hudson. Pipe Flow: A Practical and
        Comprehensive Guide. Hoboken, NJ, USA: John Wiley & Sons, Inc., 2012.
        https://doi.org/10.1002/9781118275276.
        :param flow_ratio: volumetric flow ratio array for branch flow to
        combined main flow (w3 / w1 in publication)
        :return: zeta resistance value
        """
        return -0.92 + 2.0 * (2.0 - self.C_xC - self.C_M) * flow_ratio \
            + ((2.0 * self.C_yC - 1.0) * (1.0 / self.diameter_ratio) ** 4.0
               + 2.0 * (self.C_xC - 1.0)) * flow_ratio ** 2.0

    def calc_resistance_values(self, **kwargs):
        """
        Calculates resistance values based on given correlation. Must be
        defined in concrete class implementation.
        :param kwargs:
        :return: vector including the calculated element-wise resistance values
        along the channel
        """
        # Calculate both resistance values, simpler due to using the
        # parameters self.C_xC, self.C_yC and self.C_M which could be arrays as
        # well. If only the values at the required positions should be
        # calculated, the parameters should be made function inputs as well.
        value = np.zeros(self.value.shape)
        value_dividing = \
            self.calc_dividing_junction_value(1.0 - self.main_flow_ratio)
        value_combining = \
            self.calc_combining_junction_value(1.0 - self.main_flow_ratio)
        # Determine indices for dividing or combining elements
        id_dividing = np.nonzero(self.mass_source_sign <= 0)
        id_combining = np.nonzero(self.mass_source_sign > 0)
        # Assign result values to corresponding elements
        value[id_dividing] = value_dividing[id_dividing]
        value[id_combining] = value_combining[id_combining]
        return value
