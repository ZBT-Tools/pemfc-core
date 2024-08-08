# general imports
import numpy as np
import math
from abc import ABC, abstractmethod

# local modul imports
from . import interpolation as ip, global_functions as g_func, \
    flow_resistance as fr, output_object as oo
from .fluid import fluid as fluids, evaporation_model as evap
try:
    import pemfc.src.cython.channel_heat_transfer as cht
    CHT_FOUND = True
except ModuleNotFoundError:
    CHT_FOUND = False


class Channel(oo.OutputObject1D, ABC):
    def __new__(cls, channel_dict, fluid, number=None):

        if type(fluid) is fluids.IncompressibleFluid \
                or type(fluid) is fluids.ConstantFluid:
            return super(Channel, cls).__new__(IncompressibleFluidChannel)
        elif type(fluid) in (fluids.GasMixture, fluids.CanteraGasMixture):
            return super(Channel, cls).__new__(GasMixtureChannel)
        elif type(fluid) in (fluids.TwoPhaseMixture,
                             fluids.CanteraTwoPhaseMixture):
            return super(Channel, cls).__new__(TwoPhaseMixtureChannel)
        else:
            raise NotImplementedError('Only Channel types of '
                                      'IncompressibleFluidChannel, '
                                      'GasMixtureChannel, and '
                                      'TwoPhaseMixtureChannel are '
                                      'implemented')

    def __init__(self, channel_dict, fluid, number=None):
        name = channel_dict['name']
        self.number = number
        super().__init__(name)
        self.fluid = fluid
        self._length = channel_dict.get('length', 0.1)
        self.n_nodes = len(self.fluid.density)
        self.n_ele = self.n_nodes - 1

        # element length
        self.p_out = channel_dict['p_out']
        self.pressure = self.fluid.pressure
        # self.pressure = g_func.full(self.n_nodes, self.p_out)
        # outlet and initial pressure
        temp_in = channel_dict['temp_in']
        self.temperature = self.fluid.temperature
        # self.temperature = g_func.full(self.n_nodes, self.temp_in)
        self.temp_ele = g_func.full(self.n_ele, temp_in)

        # inlet temperature
        self.tri_mtx = None
        self.id_in = None
        self.id_out = None
        self.flow_direction = channel_dict['flow_direction']

        self.pressure_recovery_factor = np.ones(self.n_ele) * 0.5

        # Geometry
        if channel_dict.keys() >= {'width', 'height'}:
            cross_shape = 'rectangular'
        elif 'diameter' in channel_dict:
            cross_shape = 'circular'
        else:
            raise KeyError('either width and height or diameter of channel '
                           'must be specified')
        self.cross_shape = channel_dict.get('cross_sectional_shape',
                                            cross_shape)
        if self.cross_shape == 'rectangular':
            self._width = channel_dict['width']
            self._height = channel_dict['height']
        elif self.cross_shape == 'circular':
            self._diameter = channel_dict['diameter']
        elif self.cross_shape == 'trapezoidal':
            self._width = channel_dict['width']
            self._height = channel_dict['height']
            self._base_width = channel_dict['base_width']
        elif self.cross_shape == 'triangular':
            self._width = channel_dict['width']
            self._height = channel_dict['height']
        else:
            raise NotImplementedError
        # Initialize geometry
        self.x, self.dx, self.dx_node = None, None, None
        self.cross_area = None
        self.surface_area = None
        self.aspect_ratio = None
        self.d_h = None
        self.base_area = None
        self.base_area_dx = None
        self.calculate_geometry()

        # Flow resistances
        self.zetas = self.create_flow_resistances(channel_dict)

        # Flow
        self.velocity = np.zeros(self.n_nodes)
        self.reynolds = np.zeros(self.n_nodes)
        self.mass_flow_total = np.zeros(self.n_nodes)
        self.volume_flow = np.zeros(self.n_nodes)
        self.g_fluid = np.zeros(self.n_nodes)

        # Heat Transfer
        self.k_coeff = np.zeros(self.n_ele)
        self.heat = np.zeros(self.n_ele)
        self.wall_temp = g_func.full(self.n_ele, temp_in)

        self.add_print_data(self.temperature, 'Fluid Temperature', 'K')
        self.add_print_data(self.wall_temp, 'Wall Temperature', 'K')
        self.add_print_data(self.pressure, 'Fluid Pressure', 'Pa')

    def create_flow_resistances(self, channel_dict):
        # Basic wall resistance
        zetas = []
        if channel_dict.get('wall_friction', True):
            zetas.append(fr.FlowResistance(self, {'type': 'WallFriction'}))

        # Resistance due to bends
        n_bends = channel_dict.get('bend_number', 0)
        zeta_bends = channel_dict.get('bend_friction_factor', 0.0)
        if n_bends > 0 and zeta_bends > 0.0:
            zeta_dict = \
                {'type': 'Constant', 'value': n_bends * zeta_bends / self.n_ele}
            zetas.append(fr.FlowResistance(self, zeta_dict))

        if 'junction_resistance_model' in channel_dict:
            zeta_dict = channel_dict['junction_resistance_model']
            zetas.append(fr.FlowResistance(self, zeta_dict))
        # Additional resistances (constant or flow splitting)
        if 'flow_resistances' in channel_dict:
            if isinstance(channel_dict['flow_resistances'], (list, tuple)):
                for zeta_dict in channel_dict['flow_resistances']:
                    zetas.append(fr.FlowResistance(self, zeta_dict))
        # if 'friction_coefficients' in channel_dict:
        #     zeta_dict = \
        #         {'type': 'Junction', 'coefficients':
        #             channel_dict['friction_coefficients']}
        #     self.zetas.append(fr.FlowResistance(self, zeta_dict))
        #
        # zeta_const = channel_dict.get('constant_friction_factor', 0.0)
        # if zeta_const > 0.0:
        #     self.zetas.append(fr.FlowResistance(self, {'type': 'Constant',
        #                                                'value': zeta_const}))
        return zetas

    def calculate_geometry(self):
        self.x = np.linspace(0.0, self._length, self.n_nodes)
        self.dx = np.diff(self.x)
        self.dx_node = np.zeros(self.n_nodes)
        self.dx_node[1:-1] = np.diff(ip.interpolate_1d(self.x))
        self.dx_node[0] = 0.5 * self.dx[0]
        self.dx_node[-1] = 0.5 * self.dx[-1]
        if self.cross_shape == 'rectangular':
            self.cross_area = self._width * self._height
            perimeter = 2.0 * (self._width + self._height)
            if self._width >= self._height:
                self.aspect_ratio = self._height / self._width
            else:
                self.aspect_ratio = self._width / self._height
        elif self.cross_shape == 'circular':
            self._width = self._diameter
            self._height = self._diameter
            self.cross_area = self._diameter ** 2.0 * 0.25 * np.pi
            perimeter = self._diameter * np.pi
            self.aspect_ratio = 1.0
        elif self.cross_shape == 'trapezoidal':
            self.cross_area = \
                (self._base_width + self._width) * self._height * 0.5
            edge = abs(self._width - self._base_width) * 0.5
            side = math.sqrt(edge ** 2.0 + self._height ** 2.0)
            perimeter = 2.0 * side + self._width + self._base_width
            if self._width >= self._height:
                self.aspect_ratio = self._height / self._width
            else:
                self.aspect_ratio = self._width / self._height
        elif self.cross_shape == 'triangular':
            self.cross_area = self._width * self._height * 0.5
            side = math.sqrt((self._width * 0.5) ** 2.0 + self._height ** 2.0)
            perimeter = 2.0 * side + self._width
        else:
            raise NotImplementedError
        self.surface_area = perimeter * self.dx
        self.d_h = 4.0 * self.cross_area / perimeter
        self.base_area = self._width * self._length
        self.base_area_dx = self._width * self.dx

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def diameter(self):
        return self.d_h

    @property
    def length(self):
        return self._length

    @width.setter
    def width(self, value):
        self._width = value
        self.calculate_geometry()

    @height.setter
    def height(self, value):
        self._height = value
        self.calculate_geometry()

    @diameter.setter
    def diameter(self, value):
        self._diameter = value
        self.calculate_geometry()

    @length.setter
    def length(self, value):
        self._length = value
        self.calculate_geometry()

    def update(self, mass_flow_in=None, mass_source=None,
               wall_temp=None, heat_flux=None, update_mass=True,
               update_flow=True, update_heat=True, update_fluid=True,
               enthalpy_source=None, **kwargs):
        if mass_flow_in is not None:
            if np.sum(mass_flow_in) < 0.0:
                id_in = int(self.id_in)
                id_out = int(self.id_out)
                self.id_in = id_out
                self.id_out = id_in
                mass_flow_in = np.abs(mass_flow_in)

        if update_mass or mass_flow_in is not None or mass_source is not None:
            self.update_mass(mass_flow_in=mass_flow_in, mass_source=mass_source)
        if update_fluid:
            self.update_fluid()
        if update_flow:
            self.update_flow()

        if update_heat:
            self.update_heat(wall_temp=wall_temp, heat_flux=heat_flux,
                             enthalpy_source=enthalpy_source,
                             channel_factor=kwargs.get('channel_factor', 1.0))

    def calculate_flow_resistance_sum(self):
        """ Update and return the sum of flow resistances (zeta-values) """
        zeta_sum = 0.0
        for zeta in self.zetas:
            zeta.update()
            zeta_sum += zeta.value
        return zeta_sum

    @abstractmethod
    def update_mass(self, mass_flow_in=None, mass_source=None,
                    update_fluid=True):
        pass

    @abstractmethod
    def update_flow(self, update_fluid=False):
        pass

    @abstractmethod
    def update_heat(self, wall_temp=None, heat_flux=None, update_fluid=True,
                    enthalpy_source=None, channel_factor=1.0):
        pass

    @abstractmethod
    def update_fluid(self):
        pass

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.volume_flow[:] = self.mass_flow_total / self.fluid.density
        self.velocity[:] = np.maximum(self.volume_flow / self.cross_area, 0.0)
        self.reynolds[:] = self.velocity * self.d_h * self.fluid.density \
            / self.fluid.viscosity

    def get_upwind_values(self, values: np.ndarray, flow_direction=None):
        if flow_direction is None:
            flow_direction = self.flow_direction
        if flow_direction > 0:
            return values[:-1]
        else:
            return values[1:]

    def get_downwind_values(self, values: np.ndarray, flow_direction=None):
        if flow_direction is None:
            flow_direction = -self.flow_direction
        else:
            flow_direction = -flow_direction
        return self.get_upwind_values(values, flow_direction)

    @property
    def flow_direction(self):
        return self._flow_direction

    @flow_direction.setter
    def flow_direction(self, flow_direction):
        if flow_direction not in (-1, 1):
            raise ValueError('Member variable flow_direction '
                             'must be either 1 or -1')
        self._flow_direction = flow_direction
        if self._flow_direction == 1:
            self.id_in = 0
            self.id_out = -1
        else:
            self.id_in = -1
            self.id_out = 0
        ones = np.zeros((self.n_ele, self.n_ele))
        ones[:] = 1.0
        if self._flow_direction == 1:
            self.tri_mtx = np.tril(ones)
        else:
            self.tri_mtx = np.triu(ones)

    @abstractmethod
    def calc_mass_balance(self, *args, **kwargs):
        pass

    def calc_heat_transfer_coeff(self):
        """
        Calculates heat transfer coefficient of channel assuming element-wise
        constant wall temperature
        (Correlations should be reviewed)
        """
        prandtl = self.fluid.viscosity * self.fluid.specific_heat / \
            self.fluid.thermal_conductivity
        d_by_l = self.d_h / self._length
        sqrt_re_pr_dbyl = np.sqrt(self.reynolds * prandtl * d_by_l)
        nu_1 = 3.66
        nu_2 = 1.66 * sqrt_re_pr_dbyl
        nu_3 = (2. / (1. + 22. * prandtl)) ** 0.166667 * sqrt_re_pr_dbyl
        nu_lam = (nu_1 ** 3. + 0.7 ** 3. + (nu_2 - 0.7) ** 3.
                  + nu_3 ** 3.) ** 0.333333

        if np.any(self.reynolds >= 2300.0):
            zeta = \
                (1.8 * np.log(self.reynolds,
                              where=self.reynolds != 0) - 1.5) ** -2.
            nu_turb = zeta / 8. * self.reynolds * prandtl \
                / (1. + 12.7 * np.sqrt(zeta / 8.)
                   * (prandtl ** 0.666667) - 1.) * (1. + d_by_l ** 0.666667)
            nusselt = \
                np.where(self.reynolds < 2300.0, nu_lam,
                         np.where(self.reynolds < 1e4,
                                  (1. - (self.reynolds - 2300.) / 7700.)
                                  * nu_lam
                                  + (self.reynolds - 2300.) / 7700. * nu_turb,
                                  nu_turb))
        else:
            nusselt = nu_lam

        ht_coeff = nusselt * self.fluid.thermal_conductivity / self.d_h
        # convection coefficient between the coolant and the channel wall
        # convection area of the channel wall
        self.k_coeff[:] = ip.interpolate_1d(ht_coeff) * self.surface_area

    def calc_heat_capacitance(self, factor=1.0):
        self.g_fluid[:] = \
            factor * self.mass_flow_total * self.fluid.specific_heat

    def calc_pressure_drop(self):
        """
        Calculates the element-wise pressure drop in the channel
        """

        if np.shape(self.velocity) != np.shape(self.fluid.density):
            raise ValueError('velocity and density arrays '
                             'must be of equal shape')
        # calculate resistance based pressure drop
        dp_zeta = np.zeros(self.n_ele)
        try:
            for zeta in self.zetas:
                zeta.update()
                dp_zeta += zeta.calc_pressure_drop()
        except FloatingPointError:
            raise FloatingPointError('check if channel geometry is '
                                     'adequate for flow conditions in {}'.
                                     format(self.name))

        # Calculate influence of dynamic pressure variation due to velocity
        # changes on static pressure drop

        # pressure_recovery_factor: k = (2 - beta) / 2,
        # default: k = 0.5 for pure Bernoulli assumption
        # as in Equation 4a in:
        # Wang, Junye. "Theory of Flow Distribution in Manifolds".
        # Chemical Engineering Journal 168, Nr. 3 (April 2011): 1331–45.
        # https://doi.org/10.1016/j.cej.2011.02.050.

        # and k = theta / 2 = (2 * beta - gamma) / 2
        # as in Equation 5 and 7 in
        # Bajura, R. A., und E. H. Jones. "Flow Distribution Manifolds".
        # Journal of Fluids Engineering 98, Nr. 4 (1. Dezember 1976): 654–65.
        # https://doi.org/10.1115/1.3448441.
        # theta = 1.05 for inlet manifold and 2.60 for outlet manifold was used
        # to compare to experimental data in above reference (equations 29 and 30)

        rho1 = self.fluid.density[:-1]
        rho2 = self.fluid.density[1:]
        v1 = self.velocity[:-1]
        v2 = self.velocity[1:]
        dp_dyn = self.pressure_recovery_factor \
            * (rho2 * v2 ** 2.0 - rho1 * v1 ** 2.0) * self.flow_direction
        return dp_zeta + dp_dyn

    def calc_pressure(self):
        """
        Calculates the static channel pressure
        """
        # density_ele = ip.interpolate_1d(self.fluid.density)
        # zeta = self.flow_resistance_sum()
        dp = self.calc_pressure_drop()
        pressure_direction = -self.flow_direction
        self.pressure[:] = self.p_out
        g_func.add_source(self.pressure, dp, pressure_direction)
        # if np.any(self.p < (self.p[self.id_out] - 1e-5)):
        #     raise ValueError('Pressure dropped below outlet pressure, please '
        #                      'check boundary conditions. Pressure calculation'
        #                      'is only adequate for incompressible flows.')

    def calc_heat_transfer(self, wall_temp=None, heat_flux=None):
        """
        Calculates heat transfer to fluid and its temperature variation
        due to the heat transfer. If wall_temp is provided, the corresponding
        heat exchange is returned. If heat_flux is provided,
        the corresponding wall_temp is returned.
        wall_temp is provided,
        :param wall_temp: 1D element-based array
        :param heat_flux: 1D element-based array
        :return: if wall_temp is provided, heat array is returned;
                 if heat_flux is provided, wall temperature is returned
        """
        g_fluid = ip.interpolate_1d(self.g_fluid)
        if wall_temp is None and heat_flux is None:
            return None
        elif wall_temp is not None and heat_flux is not None:
            raise ValueError('either wall_temp or heat_flux must be provided')
        elif wall_temp is not None:
            if np.ndim(wall_temp) == 0:
                wall_temp = g_func.full(self.temp_ele.shape, wall_temp)
            elif wall_temp.shape != self.temp_ele.shape:
                raise ValueError('wall temperature array must be element-based')
            if CHT_FOUND:
                fluid_temp, heat = \
                    cht.calc_heat_transfer(wall_temp, self.temperature,
                                           g_fluid, self.k_coeff,
                                           self.flow_direction)
            else:
                fluid_temp, heat = \
                    g_func.calc_temp_heat_transfer(wall_temp, self.temperature,
                                                   g_fluid, self.k_coeff,
                                                   self.flow_direction)
            # fluid_temp, heat = \
            #     cht.calc_heat_transfer(wall_temp, self.temperature,
            #                            g_fluid, self.k_coeff,
            #                            self.flow_direction)
            # fluid_temp, heat = \
            #     g_func.calc_temp_heat_transfer(wall_temp, self.temperature,
            #                                    g_fluid, self.k_coeff,
            #                                    self.flow_direction)

            self.wall_temp[:] = wall_temp
            self.heat[:] = heat
            self.temperature[:] = fluid_temp
            self.temp_ele[:] = ip.interpolate_1d(self.temperature)
            return heat
        elif heat_flux is not None:
            if np.ndim(heat_flux) == 0:
                heat_flux = g_func.full(self.temp_ele.shape, heat_flux)
            elif heat_flux.shape != self.temp_ele.shape:
                raise ValueError('wall temperature array must be element-based')
            heat = heat_flux * self.surface_area
            dtemp = heat / g_fluid
            self.temperature[:] = self.temperature[self.id_in]
            g_func.add_source(self.temperature, dtemp,
                              self.flow_direction, self.tri_mtx)
            self.heat[:] = heat
            self.temp_ele[:] = ip.interpolate_1d(self.temperature)
            wall_temp = self.temp_ele + heat / self.k_coeff
            self.wall_temp[:] = wall_temp
            return wall_temp
        else:
            raise ValueError

    def calc_mix_temperature(self, enthalpy_source):
        # not working correctly yet
        assert enthalpy_source.shape == self.temp_ele.shape
        if np.sum(np.abs(enthalpy_source)) > 0.0:
            if self.flow_direction == -1:
                enthalpy_in = self.g_fluid[1:] * self.temperature[1:]
                self.temperature[:-1] = 1.0 / self.g_fluid[:-1] \
                    * (enthalpy_in + self.heat + enthalpy_source)
            else:
                enthalpy_in = self.g_fluid[:-1] * self.temperature[:-1]
                self.temperature[1:] = 1.0 / self.g_fluid[:1] \
                    * (enthalpy_in + self.heat + enthalpy_source)
            self.temp_ele[:] = ip.interpolate_1d(self.temperature)


class IncompressibleFluidChannel(Channel):
    def __init__(self, channel_dict, fluid, number=None):
        super().__init__(channel_dict, fluid, number)
        self.mass_source = np.zeros(self.n_ele)

    # def update(self, mass_flow_in=None, mass_source=None,
    #            wall_temp=None, heat_flux=None,
    #            update_flow=True, update_heat=True, **kwargs):
    #     self.calc_mass_balance(mass_flow_in, mass_source)
    #     self.fluid.update(self.temp, self.p)
    #     if update_flow:
    #         self.calc_flow_velocity()
    #         self.calc_pressure()
    #     if update_heat:
    #         self.calc_heat_transfer_coeff()
    #         self.calc_heat_capacitance()
    #         self.calc_heat_transfer(wall_temp=wall_temp, heat_flux=heat_flux)

    def update_mass(self, mass_flow_in=None, mass_source=None,
                    update_fluid=True):
        self.calc_mass_balance(mass_flow_in, mass_source)
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure)

    def update_flow(self, update_fluid=False):
        self.calc_flow_velocity()
        self.calc_pressure()
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure)

    def update_heat(self, wall_temp=None, heat_flux=None, update_fluid=True,
                    enthalpy_source=None, channel_factor=1.0):
        self.calc_heat_transfer_coeff()
        self.calc_heat_capacitance(factor=channel_factor)
        self.calc_heat_transfer(wall_temp=wall_temp, heat_flux=heat_flux)
        if enthalpy_source is not None:
            self.calc_mix_temperature(enthalpy_source)
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure)

    def update_fluid(self):
        self.fluid.update(self.temperature, self.pressure)

    def calc_mass_balance(self, mass_flow_in=None, mass_source=None):
        if mass_flow_in is not None:
            self.mass_flow_total[:] = mass_flow_in
        if mass_source is not None:
            self.mass_source[:] = mass_source
        g_func.add_source(self.mass_flow_total, self.mass_source,
                          self.flow_direction, self.tri_mtx)
        self.mass_flow_total[self.mass_flow_total < 0.0] = 0.0


class GasMixtureChannel(Channel):
    def __init__(self, channel_dict, fluid, number=None):
        super().__init__(channel_dict, fluid, number)
        self.mole_flow_total = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow = np.zeros(arr_shape)
        self.mass_flow = np.zeros(arr_shape)
        self.mass_source = np.zeros((self.fluid.n_species, self.n_ele))
        self.mole_source = np.zeros((self.fluid.n_species, self.n_ele))
        self.concentration_ele = np.zeros((self.fluid.n_species, self.n_ele))

        self.add_print_data(self.mole_flow, 'Mole Flow',
                            'mol/s', self.fluid.species_names)

    # def update(self, mass_flow_in=None, mass_source=None,
    #            wall_temp=None, heat_flux=None,
    #            update_flow=True, update_heat=True, **kwargs):
    #     self.calc_mass_balance(mass_flow_in, mass_source)
    #     self.fluid.update(self.temp, self.p, self.mole_flow)
    #     if update_flow:
    #         self.calc_flow_velocity()
    #         self.calc_pressure()
    #     if update_heat:
    #         self.calc_heat_transfer_coeff()
    #         self.calc_heat_capacitance()
    #         self.calc_heat_transfer(wall_temp=wall_temp, heat_flux=heat_flux)
    #         self.fluid.update(self.temp, self.p, self.mole_flow)

    def update_mass(self, mass_flow_in=None, mass_source=None,
                    update_fluid=True):
        self.calc_mass_balance(mass_flow_in, mass_source)
        self.calc_concentration()
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure, self.mole_flow)

    def update_flow(self, update_fluid=False):
        self.calc_flow_velocity()
        self.calc_pressure()
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure, self.mole_flow)

    def update_heat(self, wall_temp=None, heat_flux=None, update_fluid=True,
                    enthalpy_source=None, channel_factor=1.0):
        self.calc_heat_transfer_coeff()
        self.calc_heat_capacitance(factor=channel_factor)
        self.calc_heat_transfer(wall_temp=wall_temp, heat_flux=heat_flux)
        if enthalpy_source is not None:
            self.calc_mix_temperature(enthalpy_source)
        # if update_fluid:
        #     self.fluid.update(self.temperature, self.pressure, self.mole_flow)

    def update_fluid(self):
        self.fluid.update(self.temperature, self.pressure, self.mole_flow)

    def calc_concentration(self, mole_fraction: np.ndarray = None):
        """
        Calculates the gas phase molar concentrations.
        """
        if mole_fraction is None:
            mole_flow = ip.interpolate_along_axis(self.mole_flow, axis=-1)
            mole_fraction = self.fluid.calc_fraction(mole_flow)
        pressure = ip.interpolate_1d(self.pressure)
        total_mol_conc = pressure \
            / (self.fluid.gas_constant * self.temp_ele)
        self.concentration_ele[:] = mole_fraction * total_mol_conc
        return self.concentration_ele

    def calc_mass_balance(self, mass_flow_in=None, mass_source=None):
        """
        Calculate mass balance in 1D channel
        :param mass_flow_in: inlet mass flow
        :param mass_source: 2D array (n_species x n_elements) of discretized
                            mass source
        :return: None
        """
        if mass_flow_in is not None:
            mass_flow_in = np.asarray(mass_flow_in)
            if mass_flow_in.shape == self.mass_flow.shape:
                # mass_flow = mass_flow_in
                mass_flow_in = mass_flow_in[:, self.id_in]
                mass_flow = \
                    g_func.fill_transposed(mass_flow_in, self.mass_flow.shape)
            elif mass_flow_in.shape == self.mass_flow_total.shape:
                mass_flow = np.outer(self.fluid.mass_fraction[:, self.id_in],
                                     mass_flow_in)
            elif mass_flow_in.shape == (self.mass_flow.shape[0],):
                mass_flow = \
                    g_func.fill_transposed(mass_flow_in, self.mass_flow.shape)
            elif np.ndim(mass_flow_in) == 0:
                mass_flow_total = np.zeros(self.mass_flow_total.shape)
                mass_flow_total[:] = mass_flow_in
                mass_flow = np.outer(self.fluid.mass_fraction[:, self.id_in],
                                     mass_flow_total)
            else:
                raise ValueError(
                    'provided mass flow cannot be converted to array'
                    ' of shape: ', self.mass_flow.shape)
            self.mass_flow[:] = mass_flow
            # self.mole_flow[:] = \
            #     (mass_flow.transpose() / self.fluid.species_mw).transpose()
        if mass_source is not None:
            if np.shape(mass_source) == (self.fluid.n_species, self.n_ele):
                self.mass_source[:] = mass_source
                self.mole_source[:] = \
                    (mass_source.transpose()
                     / self.fluid.species_mw).transpose()
            else:
                raise ValueError('shape of mass_source does not conform '
                                 'to mole_source array')
        for i in range(self.fluid.n_species):
            g_func.add_source(self.mass_flow[i], self.mass_source[i],
                              self.flow_direction, self.tri_mtx)

        self.mass_flow[self.mass_flow < 0.0] = 0.0
        self.mass_flow_total[:] = np.sum(self.mass_flow, axis=0)
        # self.mass_flow[:] = \
        #     (self.mole_flow.transpose() * self.fluid.species_mw).transpose()
        self.mole_flow[:] = \
            (self.mass_flow.transpose() / self.fluid.species_mw).transpose()
        self.mole_flow_total[:] = np.sum(self.mole_flow, axis=0)


class TwoPhaseMixtureChannel(GasMixtureChannel):
    def __init__(self, channel_dict, fluid, number=None):
        super().__init__(channel_dict, fluid, number)
        self.mole_flow_gas_total = np.zeros(self.n_nodes)
        self.mass_flow_gas_total = np.zeros(self.n_nodes)
        self.volume_flow_gas = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow_liq = np.zeros(arr_shape)
        self.mole_flow_gas = np.zeros(arr_shape)
        self.mass_flow_gas = np.zeros(arr_shape)
        self.mass_flow_liq = np.zeros(arr_shape)
        self.evaporation_rate = np.zeros(self.n_ele)
        self.evaporation_heat = np.zeros(self.n_ele)
        self.humidity_ratio = np.zeros(self.n_nodes)
        self.concentration_weights = np.ones((self.fluid.n_species, 2,
                                               self.n_ele)) * 0.5
        if 'evaporation_model' in channel_dict:
            evaporation_dict = channel_dict['evaporation_model']

        else:
            evaporation_dict = {"type": "HertzKnudsenSchrage",
                                "evaporation_coefficient": 0.37}
        self.evaporation_model = evap.EvaporationModel(evaporation_dict,
                                                       self.fluid)
        self.droplet_to_slug_transition = 0.1
        self.saturation = np.ones(self.evaporation_rate.shape) * 1e-6
        self.urf = 0.5
        self.add_print_data(self.mole_flow_gas, 'Gas Mole Flow', 'mol/s',
                            self.fluid.species_names)
        self.add_print_data(self.mole_flow_liq, 'Liquid Mole Flow',
                            'mol/s', self.fluid.species_names)
        self.add_print_data(self.saturation, 'Liquid Saturation', '-')

    def update_mass(self, mass_flow_in=None, liquid_mass_flow_in=None,
                    mass_source=None, update_fluid=True):
        self.calc_mass_balance(mass_flow_in, mass_source)
        self.calc_two_phase_flow(liquid_mass_flow_in)
        self.calc_concentration()

    def update_fluid(self):
        if np.all(self.mole_flow_gas_total > 0.0):
            mole_composition = self.mole_flow
            gas_mole_composition = self.mole_flow_gas
        else:
            mole_composition = None
            gas_mole_composition = None
        self.fluid.update(self.temperature, self.pressure,
                          mole_composition, gas_mole_composition)
        # print('test')

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.volume_flow_gas[:] = self.mass_flow_gas_total / self.fluid.density
        self.volume_flow[:] = self.volume_flow_gas
        self.velocity[:] = np.maximum(self.volume_flow_gas / self.cross_area, 0.0)
        self.reynolds[:] = self.velocity * self.d_h * self.fluid.density \
            / self.fluid.viscosity

    def calc_two_phase_flow(self, liquid_mass_flow_in=None):
        """
        Calculates the condensed phase flow and updates mole and mass fractions
        """
        id_pc = self.fluid.id_pc
        mw_pc = self.fluid.gas.species_mw[id_pc]
        if liquid_mass_flow_in is not None:
            self.mass_flow_liq[id_pc][:] = liquid_mass_flow_in
            self.mole_flow_liq[id_pc][:] = liquid_mass_flow_in / mw_pc
        # else:
        # Assuming all liquid flow is removed instantly for simplification now
        self.mass_flow_liq[id_pc][:] = 0.0
        self.mole_flow_liq[id_pc][:] = 0.0
        self.mole_flow_gas[:] = self.mole_flow
        self.mass_flow_gas[:] = self.mass_flow
        # mole_fraction_water = self.fluid.gas.mole_fraction[id_pc]
        # evaporation_rate = self.evaporation_model.calc_evaporation_rate(
        #     self.temperature, self.pressure)[0]
        # condensation_rate = -ip.interpolate_1d(evaporation_rate)
        # liquid_surface_area = self.calc_liquid_surface_area(self.saturation)
        # mass_source_liquid = condensation_rate * liquid_surface_area
        # mole_source_liquid = (mass_source_liquid / mw_liquid)
        # liquid_factor = self.fluid.humidity - 1.0
        # liquid_factor[liquid_factor < 0.0] = 0.0
        mole_flow_gas = np.copy(self.mole_flow_gas)
        mole_flow_pc_gas = np.copy(mole_flow_gas[id_pc])
        mole_fraction_pc_gas_max = (self.fluid.saturation_pressure /
                                    self.pressure * 1.0)
        mole_flow_non_pc = np.sum(np.delete(mole_flow_gas, id_pc, axis=0),
                                  axis=0)
        mole_flow_pc_gas_max = (
                mole_fraction_pc_gas_max / (1.0 - mole_fraction_pc_gas_max)
                * mole_flow_non_pc)
        concentration_pc_ratio = np.divide(
            mole_flow_pc_gas, mole_flow_pc_gas_max,
            out=np.zeros(mole_flow_pc_gas.shape),
            where=mole_flow_pc_gas_max != 0.0)
        concentration_pc_ratio_uw = self.get_upwind_values(
            concentration_pc_ratio)
        concentration_pc_ratio_dw = self.get_downwind_values(
            concentration_pc_ratio)
        indices_lower = (concentration_pc_ratio_uw < 1.0)
        indices_higher = (concentration_pc_ratio_dw > 1.0)
        indices = np.nonzero(indices_lower & indices_higher)
        concentration_averaging_weights = np.ones(
            (2, *concentration_pc_ratio_uw.shape)) * 0.5
        if indices[0].size:
            # change_lower = (1.0 - concentration_pc_ratio_uw)[indices] * 0.5
            ratio_lower = (np.ones(concentration_pc_ratio_uw.shape)[indices]
                           * 0.5) * 0.0
            ratio_higher = (concentration_pc_ratio_dw - 0.5)[indices] + 0.5
            concentration_factor_non_zero = np.asarray(
                [ratio_lower, ratio_higher])
            concentration_weights_non_zero = (
                self.fluid.gas.calc_fraction(concentration_factor_non_zero))
            concentration_averaging_weights[0][indices] = (
                concentration_weights_non_zero[0])
            concentration_averaging_weights[1][indices] = (
                concentration_weights_non_zero[1])
        self.concentration_weights[id_pc] = concentration_averaging_weights

        mole_flow_pc_gas = np.where(mole_flow_pc_gas < mole_flow_pc_gas_max,
                                    mole_flow_pc_gas, mole_flow_pc_gas_max)

        mole_flow_gas[id_pc] = mole_flow_pc_gas
        mole_fraction = self.fluid.gas.calc_fraction(mole_flow_gas)
        # humidity = (mole_fraction[id_pc] * self.pressure /
        #             self.fluid.saturation_pressure)

        self.mole_flow_gas[:] = mole_flow_gas
        self.mass_flow_gas[id_pc] = mole_flow_pc_gas * mw_pc
        # g_func.add_source(self.mass_flow_liq[id_pc], mass_source_liquid,
        #                   self.flow_direction, self.tri_mtx)
        # self.mass_flow_liq[id_pc][:] = self.mole_flow_liq[id_pc] * mw_pc

        # self.mass_flow_liq[self.mass_flow_liq < 0.0] = 0.0
        # self.mole_flow_liq[self.mole_flow_liq < 0.0] = 0.0
        self.mole_flow[:] = self.mole_flow_gas + self.mole_flow_liq
        self.mass_flow[:] = self.mass_flow_gas + self.mass_flow_liq
        self.mass_flow_gas_total[:] = np.sum(self.mass_flow_gas, axis=0)
        self.mole_flow_gas_total[:] = np.sum(self.mole_flow_gas, axis=0)
        self.mole_flow_gas[self.mole_flow_gas < 0.0] = 0.0
        self.mass_flow_gas[self.mass_flow_gas < 0.0] = 0.0

        # mole_flow_liq_total = self.mole_flow_total \
        #     * self.fluid.liquid_mole_fraction
        # mass_flow_liq_total = self.mass_flow_total \
        #     * self.fluid.liquid_mass_fraction

        # self.mole_flow_gas[:] = \
        #     self.mole_flow_gas_total * self.fluid.gas.mole_fraction
        # self.mass_flow_gas[:] = \
        #     self.mass_flow_gas_total * self.fluid.gas.mass_fraction
        # self.mole_flow_liq[:] = self.mole_flow - self.mole_flow_gas
        # self.mass_flow_liq[:] = self.mass_flow - self.mass_flow_gas
        # self.mole_flow[:] = self.mole_flow_liq + self.mole_flow_gas
        self.mole_flow_total[:] = np.sum(self.mole_flow, axis=0)
        self.mass_flow_total[:] = np.sum(self.mass_flow, axis=0)
        # self.mole_flow_gas_total[:] = \
        #     self.mole_flow_total - np.sum(self.mole_flow_liq, axis=0)
        # self.mass_flow_gas_total[:] = \
        #     self.mass_flow_total - np.sum(self.mass_flow_liq, axis=0)

        self.calc_evaporation_rate()
        self.calc_evaporation_heat()
        # self.fluid.update(self.temperature, self.pressure, self.mole_flow,
        #                   self.mole_flow_gas)
        # print('test')

    def calc_concentration(self, mole_fraction: np.ndarray = None):
        if mole_fraction is None:
            mole_flow_gas = np.asarray(
                [ip.interpolate_1d(self.mole_flow_gas[i],
                                   weights=self.concentration_weights[i])
                 for i in range(self.fluid.n_species)])

            mole_fraction = self.fluid.gas.calc_fraction(mole_flow_gas)
        concentration_ele = super().calc_concentration(mole_fraction)
        return concentration_ele

    def calc_liquid_surface_area(self, saturation):
        # diameter_max = 0.5 * self.height
        # volume_max_sphere = 4.0 / 3.0 * (diameter_max / 2.0) ** 3.0 * math.pi
        # surface_max_sphere = 4.0 * (diameter_max / 2.0) ** 2.0 * math.pi
        # volume_liquid = self.cross_area * self.dx * saturation
        # number_of_max_spheres = np.mod(volume_liquid / volume_max_sphere, 1)
        # volume_max_spheres = number_of_max_spheres * volume_max_sphere
        # volume_remaining = volume_liquid - volume_max_spheres
        # diameter_remaining = (2.0 * (3.0 * volume_remaining / (4.0 * np.pi))
        #                       ** (1.0 / 3.0))
        # surface_remaining = 4.0 * (diameter_remaining / 2.0) ** 2.0 * math.pi
        # surface_spheres = (number_of_max_spheres * surface_max_sphere
        #                    + surface_remaining)
        # surface_film = self.dx * self.height
        # surface = np.where(surface_spheres < surface_film, surface_spheres,
        #                    surface_film)
        surface_film = self.dx * self.height
        return surface_film * 10.0

    def calc_evaporation_rate(self):
        """
        Calculates the molar condensation rate of the phase change species in
        the channel.
        """
        self.evaporation_rate[:] = (self.flow_direction
                                    * np.ediff1d(self.mole_flow_gas[
                                                      self.fluid.id_pc]))
        return self.evaporation_rate

    def calc_evaporation_heat(self):
        """
        Calculates the molar condensation rate of the phase change species in
        the channel.
        """
        evaporation_rate = self.evaporation_rate
        vaporization_enthalpy = \
            self.fluid.calc_vaporization_enthalpy(self.temp_ele)
        self.evaporation_heat[:] = evaporation_rate * vaporization_enthalpy

    def calc_heat_capacitance(self, factor=1.0):
        g_fluid = \
            factor * self.mass_flow_total * self.fluid.specific_heat
        g_fluid = np.maximum(g_fluid, 1e-6)
        self.g_fluid[:] = g_fluid
