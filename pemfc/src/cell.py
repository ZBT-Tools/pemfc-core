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

        """heat conductivity along and through the cell layers"""
        # self.width_straight_channels = \
        #     self.cathode.flow_field.width_straight_channels
        self.d_area = self.cathode.gde.discretization.d_area

        # Setup membrane
        self.membrane = membrane.Membrane.create(
            membrane_dict, self.cathode.gde.discretization)
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

        # Array to link all aggregate elements indexes of 1D-solution vector
        # in groups of the corresponding layers
        self.nx = self.thermal_system.shape[0]
        self.n_layer = self.nx - 1

        # Assign layer id corresponding to material name
        self.layer_id = self.create_layer_index_dict(self.n_layer)

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
            self.nx_names.insert(-2, 'Anode BPP-BPP')

        # Current density
        self.current_density = np.zeros(
            self.electrical_system.conductance[0].shape)

        # Through-plane electrochemical (MEA) cell conductance
        self.electrochemical_conductance = np.zeros(
            self.electrical_system.conductance[0][0].shape)
        # Voltage loss over the single cell stack (bpp-to-bpp)
        self.voltage_loss = np.zeros(self.electrochemical_conductance.shape)
        self.voltage = np.zeros(self.electrochemical_conductance.shape)

        # Assign results to output data
        self.add_print_data(self.current_density[self.layer_id['membrane']],
                            'Current Density', 'A/m²')
        self.add_print_data(self.temp_layer, 'Temperature', 'K',
                            sub_names=self.nx_names[:self.nx])
        self.add_print_data(self.voltage, 'Voltage', 'K')

    @staticmethod
    def create_layer_index_dict(n_layer):
        layer_id_dict = {}
        if n_layer % 2 == 0.0:
            raise ValueError('Number of material layers in fuel cell must be '
                             'uneven')
        else:
            membrane_id = n_layer // 2
            layer_id_dict['membrane'] = membrane_id
            layer_id_dict['cathode_gde'] = membrane_id - 1
            layer_id_dict['anode_gde'] = membrane_id + 1

            layer_id_dict['cathode_bpp'] = []
            layer_id_dict['anode_bpp'] = []
            for i in range(membrane_id - 1):
                layer_id_dict['cathode_bpp'].append(membrane_id - 2 - i)
                layer_id_dict['anode_bpp'].append(membrane_id + 2 + i)
        return layer_id_dict

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

    def update(self, current_density, update_channel=False,
               current_control=True, urf=None):
        """
        This function coordinates the program sequence
        """
        if urf is None:
            urf = self.urf
        current_density = (
                (1.0 - urf) * current_density + urf * self.current_density)

        self.membrane.temp[:] = (
                0.5 * (self.temp_layer[self.layer_id['cathode_gde']] +
                       self.temp_layer[self.layer_id['anode_gde']]))
        if isinstance(self.membrane, membrane.WaterTransportMembrane):
            self.cathode.w_cross_flow[:] = self.membrane.water_flux * -1.0
            self.anode.w_cross_flow[:] = self.membrane.water_flux

        self.cathode.update(current_density[self.layer_id['cathode_gde']],
                            update_channel=update_channel,
                            current_control=current_control)
        self.anode.update(current_density[self.layer_id['anode_gde']],
                          update_channel=update_channel,
                          current_control=True)
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

            self.membrane.update(
                corrected_current_density[self.layer_id['membrane']], humidity)
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

        WARNING: For now this function is deprecated since cell voltage losses
        are calculated via the voltage difference between BPP in the
        ElectricSystem class.
        """
        self.voltage_loss[:] = \
            self.membrane.voltage_loss + self.cathode.voltage_loss + self.anode.voltage_loss
        self.voltage_alarm = self.voltage_loss >= self.e_0
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

    def calc_electrochemical_conductance(self, current_density):
        """
        Calculates the element-wise area-specific electrochemical conductance
        in through-plane direction (x-direction). This inverse resistance
        includes the charge transfer resistances of the electrode and the
        ionic resistance of the membrane at the current configuration.
        This conductance is assigned to the central layer (the membrane layer)
        of the conductance array layout assuming an uneven number of layers.
        """
        current = current_density * self.membrane.discretization.d_area
        electrochemical_resistance = (
            (self.cathode.voltage_loss
             + self.membrane.voltage_loss
             + self.anode.voltage_loss)
            / current)
        self.electrochemical_conductance[:] = (
                1.0 / electrochemical_resistance)
