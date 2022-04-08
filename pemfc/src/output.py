# general imports
import numpy as np
import os
import shutil
from itertools import cycle, islice
import matplotlib
import matplotlib.pyplot as plt
import timeit
import json
import sys

# local module imports
from . import interpolation as ip
from . import global_functions as g_func
from . import stack as stack
# from ..data import input_dicts

# configure backend here
main_name = sys.argv[0]
if 'main_app.py' in main_name:
    matplotlib.use('TkAgg')

# globals
FONT_SIZE = 14
NUMBER_SIZE = 14
MARKER_SIZE = 5.0
LINE_WIDTH = 1.0
FIG_DPI = 150
FIG_SIZE = (6.4, 4.8)


class Output:

    def __init__(self, dict_output):
        self.save_csv = dict_output['save_csv']
        # switch to save the csv data
        self.save_plot = dict_output['save_plot']
        # switch to save the plot data
        self.show_loss = dict_output['show_loss']
        # switch to show the single voltage losses in the u-i-graph
        self.delimiter = ','
        self.csv_format = '%.9e'
        # object of the class Stack
        self.output_dir = dict_output.get('directory',
                                          os.path.join(os.getcwd(), 'output'))
        self.case_name = None
        if not os.path.exists(self.output_dir):
            try:
                original_umask = os.umask(0)
                # print(self.output_dir)
                os.makedirs(self.output_dir)  # , 0o0777)
            finally:
                os.umask(original_umask)

        # self.clean_directory(self.output_dir)

        # if not os.path.exists(self.output_dir):
        #     os.makedirs(self.output_dir)
            # shutil.rmtree(self.output_dir, ignore_errors=True)

    @staticmethod
    def clean_directory(directory):
        for file in os.listdir(directory):
            file_path = os.path.join(directory, file)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path, ignore_errors=True)
            except Exception as e:
                print(e)

    @staticmethod
    def set_ax_properties(ax, **kwargs):
        fontsize = kwargs.get('fontsize', FONT_SIZE)
        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'], fontsize=fontsize)
        if 'ylabel' in kwargs:
            ax.set_ylabel(kwargs['ylabel'], fontsize=fontsize)
        if 'margins' in kwargs:
            ax.margins(x=kwargs['margins'][0], y=kwargs['margins'][1])
        if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs:
            ax.set_ylim(kwargs['ylim'])
        if 'xticks' in kwargs:
            ax.set_xticks(kwargs['xticks'])
        if 'yticks' in kwargs:
            ax.set_yticks(kwargs['yticks'])
        if 'labels' in kwargs:
            ax.legend(kwargs['labels'], fontsize=fontsize)
        if 'title' in kwargs:
            ax.set_title(kwargs['title'], fontsize=fontsize)
        return ax

    def plot_lines(self, ax, x, y, colormap=None, **kwargs):
        x = np.asarray(x)
        y = np.asarray(y)
        ny = len(y)

        if x.ndim != y.ndim:
            if x.ndim in (0, 1):
                x = np.tile(x, (ny, 1))
            else:
                raise ValueError('Outer dimension of x is not one and not '
                                 'equal to outer dimension of y')
        if y.ndim == 1:
            ax.plot(x, y, marker=kwargs.get('marker', '.'),
                    markersize=kwargs.get('markersize', MARKER_SIZE),
                    fillstyle=kwargs.get('fillstyle', 'full'),
                    linewidth=kwargs.get('linewidth', LINE_WIDTH),
                    linestyle=kwargs.get('linestyle', '-'),
                    color=kwargs.get('color', 'k'))
        else:
            if colormap is not None:
                cmap = plt.get_cmap(colormap)
                colors = cmap(np.linspace(0.0, 1.0, ny))
            else:
                colors = \
                    kwargs.get('color',
                               list(islice(cycle(['k', 'b', 'r', 'g', 'y']),
                                           ny)))
            linestyles = \
                list(islice(cycle(kwargs.get('linestyle', ['-'])), ny))
            markers = \
                list(islice(cycle(kwargs.get('marker', ['.'])), ny))
            fillstyles = \
                list(islice(cycle(kwargs.get('fillstyle', ['full'])), ny))
            for i in range(ny):
                ax.plot(x[i], y[i], marker=markers[i],
                        markersize=kwargs.get('markersize', MARKER_SIZE),
                        fillstyle=fillstyles[i],
                        linewidth=kwargs.get('linewidth', LINE_WIDTH),
                        linestyle=linestyles[i],
                        color=colors[i])
        ax.grid(True)
        ax.use_sticky_edges = False
        ax.autoscale()
        ax.set_xscale(kwargs.get('xscale', 'linear'))
        ax.set_yscale(kwargs.get('yscale', 'linear'))
        ax = self.set_ax_properties(ax, **kwargs)
        return ax

    def create_figure(self, filepath, x_array, y_array, xlabels, ylabels,
                      xlims=None, ylims=None, xticks=None, yticks=None,
                      titles=None, rows=1, cols=1, **kwargs):
        nplots = rows*cols

        def check_dims(variable, correct_single_dim=False):
            if isinstance(variable, str):
                variable = [variable]
            if not isinstance(variable, (list, tuple, np.ndarray)):
                raise TypeError('variable must be provided '
                                'as tuple, list or numpy array')
            if len(variable) != nplots:
                if correct_single_dim:
                    if nplots == 1:
                        variable = [variable]
                    else:
                        raise ValueError('variable must be sequence with '
                                         'length equivalent to number of plots')
                else:
                    raise ValueError('variable must be sequence with '
                                     'length equivalent to number of plots')
            return variable

        if rows > 2:
            figsize = kwargs.get('figsize', (FIG_SIZE[0],
                                             FIG_SIZE[1] * float(rows) / 2.0))
        else:
            figsize = kwargs.get('figsize', FIG_SIZE)
        fig = plt.figure(dpi=kwargs.get('dpi', FIG_DPI), figsize=figsize)

        x_array = np.asarray(x_array)
        y_array = check_dims(np.asarray(y_array), correct_single_dim=True)

        if len(x_array) != nplots:
            if x_array.ndim == 1:
                x_array = np.tile(x_array, (nplots, 1))
            else:
                raise ValueError('Dimension of x-array is not one and does not '
                                 'match number of plot')
        fontsize = kwargs.get('fontsize', FONT_SIZE)
        xlabels = check_dims(xlabels)
        ylabels = check_dims(ylabels)

        for i in range(nplots):
            ax = fig.add_subplot(rows, cols, i+1)
            ax = self.plot_lines(ax, x_array[i], y_array[i],
                                 xlabel=xlabels[i], ylabel=ylabels[i], **kwargs)
            if titles is not None:
                titles = check_dims(titles)
                ax.set_title(titles[i], fontsize=fontsize)
            if 'legend' in kwargs:
                legend = check_dims(kwargs['legend'], correct_single_dim=True)
                ax.legend(legend[i])
            if xlims is not None:
                xlims = check_dims(xlims, correct_single_dim=True)
                ax.set_xlim(xlims[i])
            if ylims is not None:
                xlims = check_dims(xlims, correct_single_dim=True)
                ax.set_ylim(ylims[i])
            if xticks is not None:
                xticks = check_dims(xticks, correct_single_dim=True)
                ax.set_xticks(xticks[i])
            if yticks is not None:
                xlims = check_dims(yticks, correct_single_dim=True)
                ax.set_yticks(yticks[i])
        plt.tight_layout()
        if filepath:
            fig.savefig(filepath, format=kwargs.get('fileformat', 'png'))
        return fig

    @staticmethod
    def plot(y_values, y_label, x_label, y_scale, colors,
             title, xlim_low, xlim_up, labels, path):
        if labels is not None:
            for i in range(len(y_values)):
                plt.plot(y_values[i], color=colors[i],
                         marker='.', label=labels[i])
        else:
            for i in range(len(y_values)):
                plt.plot(y_values[i], color=colors[i], marker='.')

        plt.xlabel(x_label, fontsize=FONT_SIZE)
        plt.ylabel(y_label, fontsize=FONT_SIZE)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=NUMBER_SIZE)
        plt.xlim(xlim_low, xlim_up)
        plt.tight_layout()
        plt.grid()
        if labels is not None:
            plt.legend()
        plt.savefig(os.path.join(path, title + '.png'))
        plt.close()

    @staticmethod
    def x_plot(path, x, y, x_label, y_label, x_scale='linear',
               y_scale='linear', xlim=None, ylim=None, title=None, labels=None):
        if labels is not None:
            for i in range(len(y)):
                plt.plot(x, y[i],
                         color=plt.cm.coolwarm(i / len(y)),
                         marker='.', label=labels[i])
        else:
            for i in range(len(y)):
                plt.plot(x, y[i],
                         color=plt.cm.coolwarm(i / len(y)), marker='.')

        plt.xlabel(x_label, fontsize=FONT_SIZE)
        plt.ylabel(y_label, fontsize=FONT_SIZE)
        plt.xscale(x_scale)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=NUMBER_SIZE)
        plt.autoscale(tight=True, axis='both', enable=True)
        if xlim is not None:
            plt.xlim(xlim[0], xlim[1])
        if ylim is not None:
            plt.ylim(ylim[0], ylim[1])
        plt.tight_layout()
        plt.grid()
        if labels is not None:
            plt.legend()
        plt.savefig(os.path.join(path, title + '.png'))
        plt.close()

    def write_array_to_csv(self, file_path, array, header=None,
                           separator_lines=None, mode='a'):
        with open(file_path, mode) as file:
            if header is not None:
                file.write('# ' + header + '\n')
            if separator_lines is not None:
                for i in range(len(separator_lines)):
                    a = array[i]
                    if a.ndim == 1:
                        a = a.reshape(1, a.shape[0])
                    file.write(separator_lines[i])
                    np.savetxt(file, a,
                               delimiter=self.delimiter, fmt=self.csv_format)
            else:
                np.savetxt(file, array,
                           delimiter=self.delimiter, fmt=self.csv_format)
        return file

    def write_data(self, x_values, data_array, x_label, data_name,
                   units='-', colormap='coolwarm', **kwargs):
        file_name = kwargs.get('file_name', data_name.replace(' ', '_'))
        if kwargs.get('save_plot', self.save_plot):
            y_label = data_name + ' $[' + units + ']$'
            if 'directory' in kwargs:
                directory = kwargs['directory']
            elif 'plot_dir' in kwargs:
                directory = kwargs['plot_dir']
            else:
                raise KeyError('either keyword argument directory '
                               'or plot_dir must be provided')
            file_path = os.path.join(directory, file_name + '.png')
            self.create_figure(file_path, x_values, data_array, x_label,
                               y_label, colormap=colormap, **kwargs)
        if kwargs.get('save_csv', self.save_csv):
            if 'directory' in kwargs:
                directory = kwargs['directory']
            elif 'csv_dir' in kwargs:
                directory = kwargs['csv_dir']
            else:
                raise KeyError('either keyword argument directory '
                               'or csv_dir must be provided')
            file_path = os.path.join(directory, file_name + '.csv')
            header = kwargs.pop('header', data_name + ' [' + units + ']')
            mode = kwargs.get('write_mode', 'a')
            self.write_array_to_csv(file_path, data_array,
                                    header=header, mode=mode)

    @staticmethod
    def get_data(fc_stack):
        if not isinstance(fc_stack, stack.Stack):
            raise TypeError('argument fc_stack must be of type Stack from pemfc'
                            'module')

        # def get_oo_collection_data(oo_collection, data_dict=None, **kwargs):
        #     if data_dict is None:
        #         data_dict = {}
        #     if len(oo_collection) > 1:
        #         names = kwargs.pop('names',
        #                            [item.name + ' ' + str(i + 1)
        #                             for i, item in enumerate(oo_collection)])
        #     else:
        #         names = kwargs.pop('names', [oo_collection[0].name])
        #     for i, item in enumerate(oo_collection):
        #         for p_data in item.print_data:
        #             data_dict[names[i]] = p_data
        #     return data_dict

        def get_oo_collection_data(oo_collection, data_dict=None, **kwargs):
            if data_dict is None:
                data_dict = {}
            n_items = len(oo_collection)
            names = kwargs.pop('names', [name for name
                                         in oo_collection[0].print_data[0]])
            # data_name = kwargs.pop('name', oo_collection[0].name)
            for i, (name, content) \
                    in enumerate(oo_collection[0].print_data[0].items()):
                value = content['value']
                var_array = g_func.construct_empty_stack_array(value, n_items)
                for j, item in enumerate(oo_collection):
                    var_array[j] = item.print_data[0][name]['value']
                data_dict[names[i]] = \
                    {'value': var_array,
                     'units': oo_collection[0].print_data[0][name]['units']}

            names = kwargs.pop('names', None)
            # names = kwargs.pop('names',
            #                    [[key + ' ' + sub_key for sub_key in sub_dict]
            #                     for key, sub_dict
            #                     in oo_collection[0].print_data[1].items()])
            for i, (base_name, sub_dict) \
                    in enumerate(oo_collection[0].print_data[1].items()):
                first_key = base_name if names is None else names[i][0]
                data_dict[first_key] = {}
                for j, (sub_name, content) in enumerate(sub_dict.items()):
                    value = content['value']
                    var_array = \
                        g_func.construct_empty_stack_array(value, n_items)
                    for k, item in enumerate(oo_collection):
                        var_array[k] = \
                            item.print_data[1][base_name][sub_name]['value']
                    second_key = sub_name if names is None else names[i][1]
                    data_dict[first_key][second_key] = \
                        {'value': var_array,
                         'units': oo_collection[0].
                         print_data[1][base_name][sub_name]['units']}
            return data_dict

        # Save cell values
        cells = fc_stack.cells
        xvalues = cells[0].cathode.channel.x
        xlabel = 'Channel Location'

        data = {'Channel Location':
                {'value': xvalues, 'units': 'm', 'label': xlabel},
                'Cells':
                {'value': [i + 1 for i in range(len(cells))], 'units': '-'}}
        data = get_oo_collection_data(cells, data_dict=data)
        # Save channel values
        cathode_channels = [cell.cathode.channel for cell in fc_stack.cells]
        data = get_oo_collection_data(cathode_channels, data_dict=data)
        anode_channels = [cell.anode.channel for cell in fc_stack.cells]
        data = get_oo_collection_data(anode_channels, data_dict=data)

        # Save fluid values
        cathode_fluids = [cell.cathode.channel.fluid for cell in fc_stack.cells]
        data = get_oo_collection_data(cathode_fluids, data_dict=data)
        anode_fluids = [cell.anode.channel.fluid for cell in fc_stack.cells]
        data = get_oo_collection_data(anode_fluids, data_dict=data)

        # Save fuel circuit values
        if fc_stack.n_cells > 1:
            fuel_circuits = fc_stack.fuel_circuits
            data = get_oo_collection_data(fuel_circuits, data_dict=data,
                                          names=['Cathode', 'Anode'])

        # Save coolant circuit values
        if fc_stack.coolant_circuit is not None:
            coolant_circuits = [fc_stack.coolant_circuit]
            data['Coolant Channels'] = \
                {'value': [i + 1 for i
                           in range(coolant_circuits[0].n_channels)],
                 'units': '-'}
            data = get_oo_collection_data(coolant_circuits, data_dict=data)

            cool_channels = \
                [channel for channel in coolant_circuits[0].channels]

            cool_fluids = [channel.fluid for channel in cool_channels]
            data = get_oo_collection_data(cool_fluids, data_dict=data)

        return data

    def save(self, folder_name, fc_stack):
        self.case_name = folder_name
        case_path = os.path.join(self.output_dir, folder_name)
        if not os.path.exists(case_path):
            os.makedirs(case_path)
        # else:
        #     self.clean_directory(case_path)

        if not self.save_csv and not self.save_plot:
            return None
        if not isinstance(fc_stack, stack.Stack):
            raise TypeError('argument fc_stack must be of type Stack from pemfc'
                            'module')
        csv_path = os.path.join(case_path, 'csv_data')
        plot_path = os.path.join(case_path, 'plots')
        if not os.path.exists(csv_path):
            os.makedirs(csv_path)
        # else:
        #     self.clean_directory(csv_path)
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        # else:
        #     self.clean_directory(plot_path)

        def save_oo_collection(oo_collection, x_values, x_label, **kwargs):
            # data_dict = kwargs.get('data_dict', {})
            if not hasattr(oo_collection[0], 'print_data'):
                raise TypeError
            n_items = len(oo_collection)
            for name, content in oo_collection[0].print_data[0].items():
                value = content['value']
                var_array = g_func.construct_empty_stack_array(value, n_items)
                for i, item in enumerate(oo_collection):
                    var_array[i] = item.print_data[0][name]['value']
                    # data_dict[item.name + ' ' + str(i)] =
                x = x_values
                if var_array.shape[-1] == (len(x_values) - 1):
                    x = ip.interpolate_1d(x_values)
                self.write_data(x, var_array, x_label, name, content['units'],
                                plot_dir=plot_path, csv_dir=csv_path, **kwargs)

            for base_name, sub_dict in oo_collection[0].print_data[1].items():
                for sub_name, content in sub_dict.items():
                    value = content['value']
                    var_array = \
                        g_func.construct_empty_stack_array(value, n_items)
                    for i, item in enumerate(oo_collection):
                        var_array[i] = \
                            item.print_data[1][base_name][sub_name]['value']
                    x = x_values
                    if var_array.shape[-1] == (len(x_values) - 1):
                        x = ip.interpolate_1d(x_values)
                    name = sub_name + ' ' + base_name
                    self.write_data(x, var_array, x_label, name,
                                    content['units'], plot_dir=plot_path,
                                    csv_dir=csv_path, **kwargs)

        # Save cell values
        cells = fc_stack.cells
        xvalues = cells[0].cathode.channel.x
        xlabel = 'Channel Location [m]'
        save_oo_collection(cells, xvalues, xlabel)

        # # Save half cell values
        # cathodes = [cell.cathode for cell in fc_stack.cells]
        # save_oo_collection(cathodes, xvalues, xlabel)
        # anodes = [cell.anode for cell in fc_stack.cells]
        # save_oo_collection(cathodes, xvalues, xlabel)
        data = {'channel_location':
                {'value': xvalues, 'units': 'm', 'label': xlabel}}

        # Save channel values
        cathode_channels = [cell.cathode.channel for cell in fc_stack.cells]
        save_oo_collection(cathode_channels, xvalues, xlabel)
        anode_channels = [cell.anode.channel for cell in fc_stack.cells]
        save_oo_collection(anode_channels, xvalues, xlabel)

        # Save fluid values
        cathode_fluids = [cell.cathode.channel.fluid for cell in fc_stack.cells]
        save_oo_collection(cathode_fluids, xvalues, xlabel)

        anode_fluids = [cell.anode.channel.fluid for cell in fc_stack.cells]
        save_oo_collection(anode_fluids, xvalues, xlabel)

        # Save membrane values
        membranes = [cell.membrane for cell in fc_stack.cells]
        save_oo_collection(membranes, xvalues, xlabel)

        # Save fuel circuit values
        if fc_stack.n_cells > 1:
            fuel_circuits = fc_stack.fuel_circuits
            xvalues = [i + 1 for i in range(len(cells))]
            xlabel = 'Cell'
            save_oo_collection(fuel_circuits, xvalues, xlabel,
                               legend=['Cathode', 'Anode'],
                               file_name='Fuel_Distribution')

        # Save coolant circuit values
        if fc_stack.coolant_circuit is not None:
            coolant_circuits = [fc_stack.coolant_circuit]
            xvalues = [i + 1 for i in range(coolant_circuits[0].n_channels)]
            xlabel = 'Channel'
            save_oo_collection(coolant_circuits, xvalues, xlabel,
                               file_name='Coolant_Distribution')

            cool_channels = \
                [channel for channel in coolant_circuits[0].channels]
            xvalues = cool_channels[0].x
            xlabel = 'Channel Location [m]'
            save_oo_collection(cool_channels, xvalues, xlabel)
            # Save fluid values
            cool_fluids = [channel.fluid for channel in cool_channels]
            save_oo_collection(cool_fluids, xvalues, xlabel)

    def plot_polarization_curve(self, voltage_loss,
                                cell_voltages, current_density):
        """
        Plots the polarization curve of the given
        current densities and average stack voltages.
        """
        cd_array = np.asarray(current_density) * 1.e-4
        plt.plot(cd_array, cell_voltages, marker='.', color='k',
                 label='Simulation')
        if self.show_loss is True:
            plt.plot(cd_array, voltage_loss['membrane']['average'],
                     color='b', marker='.', label='Membrane Loss')
            plt.plot(cd_array, voltage_loss['activation']['anode']['average'],
                     color='g', marker='*', label='Anode Activation Loss')
            plt.plot(cd_array,
                     voltage_loss['activation']['cathode']['average'],
                     color='g', marker='+', label='Cathode Activation Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['CL']['anode']['average'],
                     color='y', marker='*', label='Anode CL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['CL']['cathode']['average'],
                     color='y', marker='+', label='Cathode CL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['GDL']['anode']['average'],
                     color='m', marker='*', label='Anode GDL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['GDL']['cathode']['average'],
                     color='m', marker='+', label='Cathode GDL Diff Loss')
        plt.ylabel('Voltage $[V]$', fontsize=FONT_SIZE)
        plt.xlabel('Current Density $[A/cm²]$', fontsize=FONT_SIZE)
        plt.tick_params(labelsize=NUMBER_SIZE)
        plt.grid()
        plt.legend()
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.ylim(0., 1.)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'polarization_curve.png'))
        plt.close()

    def print_global_data(self, sim, data, **kwargs):
        name = kwargs.pop('name', 'summary.txt')
        summary_file_path = os.path.join(self.output_dir, self.case_name, name)
        with open(summary_file_path, 'w') as file:
            file.write('Initialization time: {0:.4f}\n'.format(
                (sim.timing['initialization'] - sim.timing['start'])))
            file.write('Simulation time: {0:.4f}\n'.format(sim.timing[
                                                               'simulation']))
            file.write('Output time: {0:.4f}\n'.format(sim.timing['output']))
            stop_time = timeit.default_timer()
            file.write('Total time: {0:.4f}\n'.format(stop_time - sim.timing[
                'start']))
            for k, v in data.items():
                file.write('{} [{}]: '.format(k, v['units'])
                           + ('{0:' + v.get('format', '.4f')
                           + '}\n').format(v['value']))

    def save_settings(self, settings, fmt='json'):
        # if settings is None:
        #     settings = input_dicts.sim_dict
        if not isinstance(settings, dict):
            raise TypeError('must provide python dict to save settings')
        else:
            file_path = os.path.join(self.output_dir, 'settings.json')
            if fmt == 'json':
                with open(file_path, 'w') as file:
                    file.write(json.dumps(settings, indent=2))

    def save_global_results(self, results, fmt='json'):
        # if settings is None:
        #     settings = input_dicts.sim_dict
        if not isinstance(results, dict):
            raise TypeError('must provide python dict to save global results')
        else:
            file_path = os.path.join(self.output_dir, self.case_name,
                                     'summary.json')
            if fmt == 'json':
                with open(file_path, 'w') as file:
                    file.write(json.dumps(results, indent=2))
