# general imports
from abc import ABC, abstractmethod
import numpy as np

# local module imports
from . import interpolation as ip
from . import constants


class FlowResistance(ABC):
    def __new__(cls, channel, zeta_dict, **kwargs):
        zeta_type = zeta_dict.get('type', 'Constant')
        if zeta_type == 'Constant':
            return super(FlowResistance, cls).\
                __new__(ConstantFlowResistance)
        elif zeta_type == 'WallFriction':
            return super(FlowResistance, cls).\
                __new__(WallFrictionFlowResistance)
        elif zeta_type == 'Junction':
            return super(FlowResistance, cls).\
                __new__(JunctionFlowResistance)
        else:
            raise NotImplementedError

    def __init__(self, channel, zeta_dict, **kwargs):
        self.channel = channel
        self.value = None

    def update(self):
        pass

    @abstractmethod
    def calc_pressure_drop(self):
        pass


class ConstantFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.value = zeta_dict['value']

    def calc_pressure_drop(self):
        dp_node = 0.5 * self.channel.fluid.density * self.value \
                  * self.channel.velocity ** 2.0
        return ip.interpolate_1d(dp_node)


class WallFrictionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.method = zeta_dict.get('method', 'Blasius')
        self.value = np.zeros(self.channel.n_nodes)

    def update(self):
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
        dp_node = 0.5 * self.channel.fluid.density * self.value \
                  * self.channel.velocity ** 2.0
        # weighting according to dx and dx_node lengths for element-wise
        # pressure drop
        dp_node_1 = dp_node[:-1]
        dp_node_2 = dp_node[1:]
        dx_node_1 = self.channel.dx_node[:-1]
        dx_node_2 = self.channel.dx_node[1:]
        dx = self.channel.dx
        # dp = (dx_node_1 * dp_node_1 + dx_node_2 * dp_node_2) / dx
        dp = (dp_node_1 + dp_node_2) * dx / (dx_node_1 + dx_node_2)
        return dp


class JunctionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.coeffs = zeta_dict['coefficients']
        self.value = np.zeros(self.channel.n_ele)
        self.branch_combined_ratio = np.zeros(self.channel.n_ele)
        self.velocity = np.zeros(self.channel.n_ele)

    def update(self):
        if hasattr(self.channel.fluid, 'species'):
            mass_source = np.sum(self.channel.mass_source, axis=0)
        else:
            mass_source = self.channel.mass_source
        mass_flow = self.channel.mass_flow_total
        if self.channel.flow_direction == 1:
            mass_flow_in = mass_flow[:-1]
            mass_flow_out = mass_flow[1:]
            velocity_in = self.channel.velocity[:-1]
            velocity_out = self.channel.velocity[1:]
        else:
            mass_flow_out = mass_flow[:-1]
            mass_flow_in = mass_flow[1:]
            velocity_out = self.channel.velocity[:-1]
            velocity_in = self.channel.velocity[1:]
        mass_source_sign = np.sign(np.sum(mass_source))
        mass_source_abs = np.abs(mass_source)
        if mass_source_sign < 0:
            # self.branch_combined_ratio[:] = mass_source_abs / mass_flow_in
            self.branch_combined_ratio[:] = \
                np.divide(mass_source_abs, mass_flow_in,
                          out=np.zeros(mass_source_abs.shape),
                          where=mass_flow_in != 0)
            self.velocity[:] = velocity_in
        else:
            # self.branch_combined_ratio[:] = mass_source_abs / mass_flow_out
            self.branch_combined_ratio[:] = \
                np.divide(mass_source_abs, mass_flow_out,
                          out=np.zeros(mass_source_abs.shape),
                          where=mass_flow_out != 0)
            self.velocity[:] = velocity_out

        self.value[:] = \
            np.polynomial.polynomial.polyval(self.branch_combined_ratio,
                                             self.coeffs)

    def calc_pressure_drop(self):
        return 0.5 * ip.interpolate_1d(self.channel.fluid.density) * \
               self.value * self.velocity ** 2.0
