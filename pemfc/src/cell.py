# General imports
import numpy as np
import copy
import math

# Local modul imports
from . import (
    half_cell as h_c, membrane as membrane, linear_system as ls)
from .output_object import OutputObject2D


class Cell(OutputObject2D):

    def __init__(self, cell_dict, membrane_dict, half_cell_dicts,
                 channels, number=None):
        name = 'Cell'  # + str(number)
        self.number = number
        super().__init__(name)
        self.cell_dict = cell_dict

        # Set number of electrodes
        self.n_electrodes = 2

        self.first_cell = cell_dict['first_cell']
        self.last_cell = cell_dict['last_cell']

        # Number of nodes/elements along the channel
        n_nodes = channels[0].n_nodes
        self.n_ele = n_nodes - 1

        # Additional resolution in width direction in regions
        # "under land" and "under channel" if True
        self.channel_land_discretization = \
            cell_dict['channel_land_discretization']

        # Additional resolution in through-plane direction (x-axis) either
        # if channel_land_discretization is True or
        # if given as keyword argument
        if self.channel_land_discretization:
            self.additional_layer = True
        else:
            self.additional_layer = cell_dict.get('additional_layer', False)
        # self.additional_layer = True

        # Underrelaxation factor
        self.urf = cell_dict['underrelaxation_factor']

        self.e_0 = cell_dict['open_circuit_voltage']
        self.e_tn = cell_dict['thermoneutral_voltage']

        self.width = self.cell_dict['width']
        self.length = self.cell_dict['length']
        self.active_area = self.width * self.length

        half_cell_dicts = copy.deepcopy(half_cell_dicts)
        # Create half cell objects
        for i in range(len(half_cell_dicts)):
            # name = self.name + ': ' + half_cell_dicts[i]['name']
            name = half_cell_dicts[i]['name']
            half_cell_dicts[i]['name'] = name
        self.half_cells = [h_c.HalfCell(half_cell_dicts[i], cell_dict,
                                        channels[i], number=self.number)
                           for i in range(len(half_cell_dicts))]
        self.cathode = self.half_cells[0]
        self.anode = self.half_cells[1]

        # Append half cell names to output data
        for half_cell in self.half_cells:
            half_cell.channel.extend_data_names(half_cell.channel.name)
            half_cell.channel.fluid.extend_data_names(
                half_cell.channel.fluid.name)

        # self.width_straight_channels = \
        #     self.cathode.flow_field.width_straight_channels
        self.discretization = self.cathode.gde.discretization
        self.d_area = self.discretization.d_area

        # Setup membrane
        self.membrane = membrane.Membrane.create(
            membrane_dict, self.discretization)
        # memb = membrane.Constant(
        #     membrane_dict, self.cathode.gde.discretization)
        # Overall cell thickness (cell pitch)
        self.thickness = self.cathode.thickness + self.membrane.thickness \
            + self.anode.thickness

        # Cell coordinates in z-direction (stacking/current direction)
        # will be initialized correctly through stack class
        self.coords = [0.0, 0.0]

        # Layer stack
        self.layers = [
            self.cathode.bpp,
            self.cathode.gde,
            self.membrane,
            self.anode.gde,
            self.anode.bpp
        ]
        # Layer thickness stack
        self.th_layer = [layer.thickness for layer in self.layers]

        # Setup linear systems for different transport physics within cell
        self.transport_types = ['thermal', 'electrical']
        # Initializing temperatures with average channel fluid temperature
        temp_init = np.average([hc.channel.fluid.temperature
                                for hc in self.half_cells])
        init_values = [temp_init, 0.0]
        self.linear_systems = {
            name: ls.CellLinearSystem(self, name, init_values[i])
            for i, name in enumerate(self.transport_types)}
        self.thermal_system = self.linear_systems['thermal']
        self.electrical_system = self.linear_systems['electrical']

        # Set thermal boundary conditions at end plates
        if self.first_cell:
            self.thermal_system.set_layer_boundary_conditions(layer_id=0)
        if self.last_cell:
            self.thermal_system.set_layer_boundary_conditions(layer_id=-1)

        # Boolean alarm values
        self.voltage_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        # Cell thickness
        self.thickness = np.sum([layer.thickness for layer in self.layers])
        # Cell temperatures
        self.temp_layer = self.thermal_system.solution_array
        # Cell voltage
        self.voltage_layer = self.electrical_system.solution_array

        # Array to link all aggregate elements indexes of 1D-solution vector
        # in groups of the corresponding layers
        self.nx = self.thermal_system.shape[0]
        self.n_layer = self.nx - 1
        # Interface names according to temperature array
        self.nx_names = [
            'Cathode BC-BPP',
            'Cathode BPP-GDE',
            'Cathode GDE-MEM',
            'Anode MEM-GDE',
            'Anode GDE-BPP',
            'Anode BPP-BC']
        if self.additional_layer:
            self.nx_names.insert(1, 'Cathode BPP-BPP')
            self.nx_names.insert(-1, 'Anode BPP-BPP')

        # Assign layer id corresponding to material name
        self.layer_id, self.interface_id = self.create_layer_index_dict()

        # Current density
        self.current_density = np.zeros(
            self.electrical_system.conductance[0].shape)

        # Through-plane electrochemical (MEA) cell conductance
        self.electrochemical_conductance = np.zeros(
            self.electrical_system.conductance[0][0].shape)
        # Voltage loss over the single cell stack (bpp-to-bpp)
        self.voltage_loss = np.zeros(self.electrochemical_conductance.shape)
        self.voltage = np.zeros(self.electrochemical_conductance.shape)
        # Heat source array
        self.explicit_heat_sources = np.zeros(self.temp_layer.shape)

        # Assign results to output data
        self.add_print_data(self.current_density[self.layer_id['membrane']],
                            'Current Density', 'A/m²')
        self.add_print_data(self.temp_layer, 'Temperature', 'K',
                            sub_names=self.nx_names[:self.nx])
        self.add_print_data(self.voltage, 'Voltage', 'V')

    def create_layer_index_dict(self):
        if self.n_layer % 2 == 0.0:
            raise ValueError('Number of material layers in fuel cell must be '
                             'uneven')
        else:
            interface_keys = [name.replace('-', '_').replace(' ', '_').lower()
                              for name in self.nx_names]
            interface_id_dict = {key: value for value, key in enumerate(
                interface_keys)}

            layer_id_dict = {}
            membrane_id = self.n_layer // 2
            layer_id_dict['membrane'] = membrane_id
            layer_id_dict['cathode_gde'] = membrane_id - 1
            layer_id_dict['anode_gde'] = membrane_id + 1

            layer_id_dict['cathode_bpp'] = []
            layer_id_dict['anode_bpp'] = []
            for i in range(membrane_id - 1):
                layer_id_dict['cathode_bpp'].append(membrane_id - 2 - i)
                layer_id_dict['anode_bpp'].append(membrane_id + 2 + i)
        return layer_id_dict, interface_id_dict

    def calc_ambient_conductance(self, alpha_amb: float,
                                 th_layer_amb: list[np.ndarray]):
        """
        :param alpha_amb: heat transfer coefficient for free or forced
        convection to the ambient
        :return: discretized conductance for each layer of cell based on
        external surface
        """
        dy = self.cathode.discretization.dx[0].flatten(order='F')
        k_amb = np.outer(th_layer_amb, dy) \
            * alpha_amb * self.cathode.flow_field.external_surface_factor
        return k_amb

    def update(self, current_density, update_channel=True,
               current_control=True, urf=None):
        """
        This function coordinates the program sequence
        """
        if urf is None:
            urf = self.urf
        current_density = (
                (1.0 - urf) * current_density + urf * self.current_density)

        self.membrane.temp[:] = (
                0.5 * (self.temp_layer[self.interface_id['cathode_gde_mem']] +
                       self.temp_layer[self.interface_id['anode_mem_gde']]))
        if isinstance(self.membrane, membrane.WaterTransportMembrane):
            self.cathode.w_cross_flow[:] = self.membrane.water_flux * -1.0
            self.anode.w_cross_flow[:] = self.membrane.water_flux

        # Update electrochemical heat sources
        self.calc_electrochemical_heat_sources(
            current_density[self.layer_id['membrane']])
        # Calculate heat flux at MEA-GDL interfaces
        heat_flux_cat_gdl, heat_flux_ano_gdl = self.calc_heat_flux()

        mea_current_density_prev = (
            self.current_density[self.layer_id['membrane']])
        mea_current_density = current_density[self.layer_id['membrane']]
        cathode_temperature = np.stack(
            [self.temp_layer[self.interface_id['cathode_bpp_gde']],
             self.temp_layer[self.interface_id['cathode_gde_mem']]], axis=0)
        self.cathode.update(mea_current_density,
                            cathode_temperature,
                            update_channel=update_channel,
                            current_control=current_control,
                            heat_flux=heat_flux_cat_gdl)
        anode_temperature = np.stack(
            [self.temp_layer[self.interface_id['anode_mem_gde']],
             self.temp_layer[self.interface_id['anode_gde_bpp']]], axis=0)
        self.anode.update(mea_current_density,
                          anode_temperature,
                          update_channel=update_channel,
                          current_control=True,
                          heat_flux=heat_flux_ano_gdl)

        # Correct common membrane water cross-flow in case of one-side channel
        # dry-out
        #TODO: Update water cross flow, however only calculate differences to
        # previous water cross flow (mole/m²-s), where mole_sources (mole/s)
        # are smaller than previous water cross flow calculations (all in
        # magnitude)

        if self.cathode.electrochemistry.corrected_current_density is not None:
            corrected_current_density = \
                self.cathode.electrochemistry.corrected_current_density
        else:
            corrected_current_density = current_density
        if self.anode.break_program or self.cathode.break_program:
            self.break_program = True
        else:
            humidity = np.asarray([self.cathode.calc_humidity(),
                                   self.anode.calc_humidity()])
            liquid_pressure = np.asarray([self.cathode.calc_liquid_pressure(),
                                          self.anode.calc_liquid_pressure()])

            self.membrane.update(
                corrected_current_density[self.layer_id['membrane']], humidity,
                liquid_pressure=liquid_pressure)
            # self.calc_voltage_loss()
            self.calc_electrochemical_conductance(
                corrected_current_density[self.layer_id['membrane']])
            self.current_density[:] = current_density
            self.voltage[:] = self.e_0 - self.voltage_loss

    def calc_voltage_loss(self):
        """
        Calculates the cell voltage loss consisting of the half cell and
        membrane voltage losses.If the cell voltage loss is larger
        than the open circuit cell voltage, the cell voltage is set to zero.

        WARNING: For now this function is only used for electrochemical
        resistance calculation. Global cell voltage losses
        are calculated via the voltage difference between BPP in the
        ElectricSystem class.
        """
        voltage_loss = (self.membrane.voltage_loss
                        + self.cathode.voltage_loss
                        + self.anode.voltage_loss)
        self.voltage_alarm = voltage_loss >= self.e_0
        return voltage_loss
        # self.v_loss[:] = np.minimum(self.v_loss, self.e_0)

    def update_voltage_loss(self, v_loss):
        """
        Function to scale internally calculated cell voltage loss with
        externally calculated overall cell voltage loss.
        Args:
            v_loss: numpy array of externally calculated local
            voltage losses
        Returns: None
        """
        pass
        # correction_factor = v_loss / self.v_loss
        # self.v_loss[:] *= correction_factor
        # for electrode in self.half_cells:
        #     electrode.v_loss[:] *= correction_factor
        # self.membrane.v_loss[:] *= correction_factor

    def calc_electrochemical_heat_sources(self, current_density: np.ndarray):
        current = current_density * self.d_area
        half_ohmic_heat_membrane = (
                0.5 * self.membrane.omega * np.square(current))

        # Cathode gde-mem source
        cathode_source = np.copy(half_ohmic_heat_membrane)
        idx = self.interface_id['cathode_gde_mem']
        v_loss = np.minimum(self.e_0, self.cathode.voltage_loss)
        v_loss[v_loss < 0.0] = 0.0
        cathode_reaction_heat = (self.e_tn - self.e_0 + v_loss) * current
        cathode_source += cathode_reaction_heat
        self.explicit_heat_sources[idx] = cathode_source

        # Anode gde-mem source
        idx = self.interface_id['anode_mem_gde']
        anode_source = np.copy(half_ohmic_heat_membrane)
        v_loss = np.minimum(self.e_0, self.anode.voltage_loss)
        v_loss[v_loss < 0.0] = 0.0
        anode_reaction_heat = v_loss * current
        anode_source += anode_reaction_heat
        self.explicit_heat_sources[idx] = anode_source

    def calc_heat_flux(self):
        # Retrieve layer and interface ids
        idx_mem = self.layer_id['membrane']
        idx_cat_gde_mem = self.interface_id['cathode_gde_mem']
        idx_ano_mem_gde = self.interface_id['anode_mem_gde']

        # Conductance [W/K] of membrane in x-direction (index: 0)
        cond_mem = self.thermal_system.conductance[0][idx_mem]

        # Temperatures [K] at cathode and anode gde-membrane interfaces
        temp_cat_gde_mem = self.temp_layer[idx_cat_gde_mem]
        temp_ano_mem_gde = self.temp_layer[idx_ano_mem_gde]

        # Heat sources
        heat_source_cat = self.explicit_heat_sources[idx_cat_gde_mem]
        heat_source_ano = self.explicit_heat_sources[idx_ano_mem_gde]

        # Heat [W] in membrane from cathode to anode side
        heat_membrane = - cond_mem * (temp_ano_mem_gde - temp_cat_gde_mem)

        # Energy balance at mem-gde interfaces to calculate heat flux [W/m²]
        heat_flux_cat_gdl = ((heat_source_cat - heat_membrane) /
                             self.cathode.discretization.d_area)
        heat_flux_ano_gdl = ((heat_source_ano + heat_membrane) /
                             self.anode.discretization.d_area)
        return heat_flux_cat_gdl, heat_flux_ano_gdl

    def calc_electrochemical_conductance(self, current_density: np.ndarray):
        """
        Calculates the element-wise area-specific electrochemical conductance
        in through-plane direction (x-direction). This inverse resistance
        includes the charge transfer resistances of the electrode and the
        ionic resistance of the membrane at the current configuration.
        This conductance is assigned to the central layer (the membrane layer)
        of the conductance array layout assuming an uneven number of layers.
        """
        current = current_density * self.membrane.discretization.d_area
        mea_voltage_loss = self.calc_voltage_loss()
        electrochemical_resistance = mea_voltage_loss / current
        self.electrochemical_conductance[:] = (
                1.0 / electrochemical_resistance)
