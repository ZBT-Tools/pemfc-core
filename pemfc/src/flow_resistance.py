# general imports
from abc import ABC
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
        self.zeta_const = zeta_dict.get('value', 0.0)
        self.factor = zeta_dict['factor']
        self.value = np.zeros(self.channel.n_ele)

    def update(self):
        pass
