import numpy as np
from . import output_object as oo


class SolidLayer(oo.OutputObject):
    def __init__(self, layer_dict, dx):
        self.name = layer_dict.get('name', 'unnamed')
        super().__init__(self.name)
        self.dx = dx
        self.layer_dict = layer_dict
        self.thickness = layer_dict['thickness']
        self.width = layer_dict['width']
        self.length = layer_dict['length']
        self.porosity = layer_dict.get('porosity', 0.0)
        self.bruggeman_exponent = layer_dict.get('bruggeman_exponent', 1.5)
        self.thermal_conductivity = \
            layer_dict.get('thermal_conductivity', 0.0)
        self.electrical_conductivity = \
            layer_dict.get('electrical_conductivity', 0.0)

        self.area = self.width * self.length
        self.area_dx = self.width * self.dx
        self.thermal_conductance = \
            self.calc_conductance(self.thermal_conductivity)
        self.electrical_conductance = \
            self.calc_conductance(self.electrical_conductivity)

    def calc_conductance(self, conductivity, effective=True):
        if np.ndim(conductivity) == 0:
            conductance_x = \
                self.width * self.thickness * conductivity / self.dx
            conductance_z = self.area_dx * conductivity / self.thickness
        elif len(conductivity) == 2:
            conductance_x = \
                self.width * self.thickness * conductivity[1] / self.dx
            conductance_z = self.area_dx * conductivity[0] / self.thickness
        else:
            raise ValueError('conductivity must be either single scalar or '
                             'an iterable with two entries')
        if effective:
            conductance_z *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_x *= (1.0 - self.porosity) ** self.bruggeman_exponent
        return np.asarray([conductance_z, conductance_x])

    def calc_voltage_loss(self, current_density, area=None, **kwargs):
        if area is None:
            current = current_density * self.area_dx
        else:
            current = current_density * area
        return current / self.electrical_conductance[0]


