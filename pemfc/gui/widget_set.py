# global imports
import tkinter as tk
from tkinter import ttk
from abc import ABC, abstractmethod
import numpy as np

# local imports
from pemfc.src import global_functions as gf
from . import base
from . import entry_value
from . import button


class WidgetSetFactory:

    def create(self, frame, **kwargs):
        widget_type = kwargs.pop('type', None)
        if widget_type == 'EntrySet':
            return self.create_entry_set(frame, **kwargs)
        elif widget_type == 'CheckButtonSet':
            return MultiCheckButtonSet(frame, **kwargs)
        elif widget_type == 'Label':
            return Label(frame, **kwargs)
        elif widget_type == 'OptionMenuSet':
            return OptionMenuSet(frame, **kwargs)
        elif widget_type == 'ComboboxSet':
            return ComboboxSet(frame, **kwargs)
        elif widget_type == 'EntryButtonSet':
            return EntryButtonSet(frame, **kwargs)
        else:
            raise NotImplementedError('type of WidgetSet not implemented')

    @staticmethod
    def create_entry_set(frame, **kwargs):
        dimensions = kwargs.get('dimensions', None)
        if dimensions is None:
            return MultiEntrySet(frame, **kwargs)
        else:
            return DimensionedEntrySet(frame, **kwargs)


class Label(base.Base):

    def __init__(self, frame, label, **kwargs):

        name = label.lower().strip(':')
        super().__init__(frame, name, **kwargs)
        self.frame = frame
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)

        kwargs['text'] = label
        self.label = tk.Label(frame, **kwargs)
        # self.label.grid(row=self.row, column=self.column, padx=self.padx,
        #                 pady=self.pady, sticky=kwargs.pop('sticky', 'W'))

    def set_grid(self, widget=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        if widget is None:
            widget = self.label
        self._set_grid(widget, row=row, column=column, **kwargs)
        return row, column

    def _get_values(self, get_object=False):
        values = {'gui_name': self.name}
        if self.sim_name is not None:
            values['sim_name'] = self.sim_name
        if get_object:
            values['object'] = self
        return values

    # @abstractmethod
    def get_values(self, get_object=False):
        return self._get_values(get_object=get_object)


class MultiWidgetSet(Label, ABC):

    WIDTH = 10

    def __init__(self, frame, label, number=None, value=None, **kwargs):
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', None)
        self.set_sticky(**kwargs)
        self.entry_value_factory = entry_value.EntryValueFactory()
        self.widgets = []
        self.number, self.init_value = self.get_number(number, value)
        self.columnspan = gf.ensure_list(self.columnspan)
        n_columnspan = self.number + 1
        if len(self.columnspan) == 1:
            self.columnspan = \
                [self.columnspan[0] for i in range(n_columnspan)]
        elif len(self.columnspan) < n_columnspan:
            columnspan_end = \
                [1 for i in range(n_columnspan - len(self.columnspan))]
            self.columnspan.extend(columnspan_end)
        self.start_columns = []
        # self.shape = None

    @staticmethod
    def get_number(number, value, dtype=None):
        if value is not None:
            # number = len(value)
            length = 1 if number is None else number
            value = gf.ensure_list(value, length=length)
            value = np.asarray(value, dtype=dtype).flatten()
            if number is None:
                number = len(value)
            # value = gf.ensure_list(value, length=number)
            # self.shape = gf.dim(value)

            return number, value
        else:
            return number, value

    @abstractmethod
    def create_widgets(self, frame, number, value, **kwargs):
        pass

    def set_grid(self, widgets=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        kwargs.pop('sticky', None)
        sticky = gf.ensure_list(self.sticky)
        row, column = super().set_grid(row=row, column=column,
                                       columnspan=self.columnspan[0],
                                       sticky=sticky[0], **kwargs)
        if widgets is None:
            widgets = self.widgets
        for i, widget in enumerate(widgets):
            column += self.columnspan[i]
            super().set_grid(widget=widget, row=row, column=column,
                             columnspan=self.columnspan[i + 1],
                             sticky=sticky[-1], **kwargs)
            # self.frame.rowconfigure(row, weight=1)
        return row, column

    def get_tk_values(self, tk_objects, get_object=False):
        values = super().get_values(get_object=get_object)
        if len(tk_objects) > 1:
            values['value'] = []
            for item in tk_objects:
                values['value'].append(
                    self.entry_value_factory.create(item.get(),
                                                    self.dtype, self))
        else:
            values['value'] = \
                self.entry_value_factory.create(tk_objects[0].get(),
                                                self.dtype, self)
        return values

    def get_values(self, get_object=False):
        return self.get_tk_values(self.widgets, get_object=get_object)

    def set_sticky(self, **kwargs):
        sticky = kwargs.pop('sticky', ['W', 'E'])
        if not isinstance(sticky, (list, tuple)):
            sticky = [sticky, 'E']
        self.sticky = sticky

    def set_tk_values(self, tk_objects, value, index=0):
        if isinstance(value, (list, tuple, np.ndarray)):
            for i, widget in enumerate(tk_objects):
                if i > len(value) - 1:
                    break
                state = widget.cget('state')
                widget.config(state='normal')
                widget.delete(0, tk.END)
                widget.insert(0,
                              self.entry_value_factory.create(value[i],
                                                              self.dtype,
                                                              self).value)
                widget.config(state=state)
        else:
            widget = tk_objects[index]
            state = widget.cget('state')
            widget.config(state='normal')
            widget.delete(0, tk.END)
            widget.insert(0,
                          self.entry_value_factory.create(value,
                                                          self.dtype,
                                                          self).value)
            widget.config(state=state)

    def set_values(self, values, index=0):
        self.set_tk_values(self.widgets, values, index=index)


class MultiEntrySet(MultiWidgetSet):

    def __init__(self, frame, label, number=None, value=None, **kwargs):
        justify = kwargs.pop('justify', 'right')
        width = kwargs.pop('width', self.WIDTH)
        super().__init__(frame, label, number=number, value=value, **kwargs)
        number, value = self.get_number(number, value)
        # dtype = 'float'
        # if value is not None:
        #     dtype = str(type(value[0]))
        self.dtype = kwargs.pop('dtype', 'float')
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)
        self.shape = (number, 0)
        kwargs['width'] = width
        kwargs['justify'] = justify
        self.create_widgets(frame, number, value, **kwargs)

    def create_widgets(self, frame, number, value, **kwargs):
        for i in range(number):
            entry = tk.Entry(frame, **kwargs)
            # entry.grid(row=self.row, column=self.column + 1 + i,
            #            padx=self.padx, pady=self.pady)
            entry.delete(0, -1)
            if value is not None:
                entry.insert(0, value[i])
            self.widgets.append(entry)


class DimensionedEntrySet(MultiEntrySet):
    def __init__(self, frame, label, number=None, dimensions='-',
                 value=None, **kwargs):
        super().__init__(frame, label, number=number, value=value, **kwargs)
        kwargs['text'] = dimensions
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)

        self.dimensions = tk.Label(frame, **kwargs)
        # self.dimensions.grid(row=self.row, column=self.column + number + 1,
        #                      padx=kwargs.get('padx', self.PADX),
        #                      pady=kwargs.get('pady', self.PADY))
        if len(self.columnspan) < (self.number + 2):
            self.columnspan.append(1)

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        row, column = super().set_grid(row=row, column=column, **kwargs)
        column += self.columnspan[self.number]
        self._set_grid(self.dimensions, row=row, column=column,
                       columnspan=self.columnspan[self.number + 1],
                       sticky='W', **kwargs)
        return row, column


class MultiCommandWidgetSet(MultiWidgetSet, ABC):

    def __init__(self, frame, label, **kwargs):
        super().__init__(frame, label, **kwargs)
        self.command_list = None
        self.execute = True

    def get_commands(self, command, number):
        self.command_list = [None for i in range(number)]
        if command is not None:
            commands = self.set_commands(command)
            len_commands = len(commands)
            self.command_list[:len_commands] = commands

    @abstractmethod
    def set_commands(self, command):
        pass

    def call_commands(self):
        if self.execute:
            if isinstance(self.command_list, list):
                for command in self.command_list:
                    if callable(command):
                        command()

    @staticmethod
    def call_object_method(obj, func, **kwargs):
        if isinstance(kwargs, dict):
            getattr(obj, func)(**kwargs)
        else:
            getattr(obj, func)()

    def call_widgets_methods(self, item_list, func, kwargs=None):
        for item in item_list:
            widget = self.frame.widget_grid[item[0]][item[1]]
            # self.call_commands()
            if isinstance(widget, tk.Widget):
                # self.call_object_method(widget, func, **kwargs)
                if isinstance(kwargs, dict):
                    getattr(widget, func)(**kwargs)
                    # self.call_object_method(widget, func, **kwargs)
                else:
                    getattr(widget, func)()
                    # self.call_object_method(widget, func)


class MultiCheckButtonSet(MultiCommandWidgetSet):
    def __init__(self, frame, label, number=None, value=None, **kwargs):
        command = kwargs.pop('command', None)

        super().__init__(frame, label, number=number, value=value, **kwargs)

        self.dtype = kwargs.pop('dtype', 'boolean')
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)
        self.check_vars = None
        # if value is not None:
        #     value = gf.ensure_list(value, length=number)
        number, value = self.get_number(number, value)
        self.get_commands(command, number)
        self.create_widgets(frame, number, value, **kwargs)

    def create_widgets(self, frame, number, value, **kwargs):
        self.check_vars = []
        for i in range(number):
            check_var = tk.BooleanVar()
            self.check_vars.append(check_var)
            check_button = \
                tk.Checkbutton(frame, variable=check_var, onvalue=True,
                               offvalue=False, command=self.command_list[i],
                               takefocus=0, **kwargs)
            if value is not None and value[i]:
                check_button.select()
                # check_button.invoke()
            self.widgets.append(check_button)

    def set_commands(self, commands_dict):
        function = commands_dict.pop('function', None)
        if function == 'set_visibility':
            arg_list = commands_dict['args']
            command_list = []
            for i, args in enumerate(arg_list):
                command_list.append(lambda arg1=i, arg2=args:
                                    self.set_visibility(arg1, arg2))
            commands_dict['function'] = self.set_visibility
            return command_list
        elif function == 'set_status':
            arg_list = commands_dict['args']
            command_list = []
            for i, args in enumerate(arg_list):
                command_list.append(lambda arg1=i, arg2=args:
                                    self.set_status(arg1, arg2))
            commands_dict['function'] = self.set_status
            return command_list
        else:
            print(function)
            raise NotImplementedError

    def widget_connector(self, widget_id, grid_list, func1, func2=None,
                         kwargs1=None, kwargs2=None):
        check_var = self.check_vars[widget_id].get()
        if check_var:
            # for item in grid_list:
            #     widget = self.frame.widget_grid[item[0]][item[1]]
            #     if isinstance(widget, tk.Widget):
            #         self.call_object_method(widget, func1, **kwargs1)
            self.call_widgets_methods(grid_list, func1, kwargs1)
        # if not self.frame.initialize:
            for item in grid_list:
                widget = self.frame.widget_grid[item[0]][item[1]]
                # self.call_commands()
                if isinstance(widget, ttk.Combobox):
                    # widget.set(widget.get())
                    widget.event_generate('<<ComboboxSelected>>')
        else:
            if func2 is None:
                func2 = func1
            # for item in grid_list:
            #     widget = self.frame.widget_grid[item[0]][item[1]]
            #     if isinstance(widget, tk.Widget):
            #         self.call_object_method(widget, func2, **kwargs2)
            self.call_widgets_methods(grid_list, func2, kwargs2)

    def set_visibility(self, widget_id, grid_list):
        self.widget_connector(widget_id, grid_list, 'grid', 'grid_remove')

    def set_status(self, widget_id, grid_list):
        self.widget_connector(widget_id, grid_list, 'config',
                              kwargs1={'state': 'normal'},
                              kwargs2={'state': 'disable'})

    def get_values(self, get_object=False):
        return super().get_tk_values(self.check_vars, get_object=get_object)

    def set_grid(self, **kwargs):
        sticky = kwargs.pop('sticky', 'WE')
        super().set_grid(sticky=sticky, **kwargs)

    def set_values(self, values, index=0):
        if isinstance(values, (list, tuple, np.ndarray)):
            for i, chkvar in enumerate(self.check_vars):
                if i > len(values) - 1 or i > len(self.check_vars) - 1:
                    break
                value = \
                    self.entry_value_factory.create(values[i],
                                                    self.dtype, self).value
                if chkvar.get() != value:
                    self.widgets[i].invoke()
                # chkvar.set(value)
                # self.widgets[i].invoke()
        else:
            chkvar = self.check_vars[index]
            value = \
                self.entry_value_factory.create(values, self.dtype, self).value
            if chkvar.get() != value:
                self.widgets[index].invoke()
            # chkvar.set(value)
            # self.widgets[index].invoke()


class OptionMenuSet(MultiCommandWidgetSet):
    def __init__(self, frame, label, number=None, **kwargs):
        value = kwargs.pop('options', [])
        super().__init__(frame, label, number=number, value=value, **kwargs)
        self.dtype = kwargs.pop('dtype', 'string')
        self.option_vars = []
        number, value = self.get_number(number, value)
        self.create_widgets(frame, number, value, **kwargs)

    def create_widgets(self, frame, number, value, **kwargs):
        for i in range(number):
            option_var = tk.StringVar()
            self.option_vars.append(option_var)
            option_menu = tk.OptionMenu(self.frame, option_var, *value)
            self.widgets.append(option_menu)

    def set_commands(self, command):
        pass

    def set_values(self, values, index=0):
        if isinstance(values, (list, tuple, np.ndarray)):
            for i in range(len(self.widgets)):
                if i > len(values) - 1 or i > len(self.widgets) - 1:
                    break
                option_list = []
                menu = self.widgets[i]['menu']
                for j in range(menu.index("end") + 1):
                    option_list.append(menu.entrycget(j, "label"))
                if values[i] in option_list:
                    self.option_vars[i].set(values[i])
                else:
                    raise ValueError('value not found in OptionMenu')
        else:
            option_list = []
            menu = self.widgets[index]['menu']
            for i in range(menu.index("end") + 1):
                option_list.append(menu.entrycget(i, "label"))
            if values[index] in option_list:
                self.option_vars[index].set(values)
            else:
                raise ValueError('value not found in OptionMenu')


class ComboboxSet(MultiCommandWidgetSet):
    def __init__(self, frame, label, number=None, value=None, **kwargs):
        # value = kwargs.pop('options', [])
        command = kwargs.pop('command', None)
        justify = kwargs.pop('justify', 'left')
        width = kwargs.pop('width', self.WIDTH)
        super().__init__(frame, label, number=number, value=value, **kwargs)
        self.values = []
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)
        self.dtype = kwargs.pop('dtype', 'string')
        # number, value = self.get_number(number, value, dtype=object)
        self.get_commands(command, number)
        kwargs['width'] = width
        kwargs['justify'] = justify
        self.create_widgets(frame, number, value, **kwargs)

    def create_widgets(self, frame, number, value, **kwargs):
        width = kwargs.pop('width', self.WIDTH)
        for i in range(number):
            self.values.append(value)
            combobox = ttk.Combobox(self.frame, values=value, width=width,
                                    **kwargs)
            combobox.current(0)
            combobox.bind('<<ComboboxSelected>>', self.command_list[i])
            self.widgets.append(combobox)

    def set_commands(self, commands_dict):
        function = commands_dict.pop('function', None)
        if function == 'show_connected_widgets':
            arg_list = commands_dict['args']
            command_list = []
            for i, args in enumerate(arg_list):
                command_list.append(lambda arg1=i, arg2=args:
                                    self.show_connected_widgets(arg1, arg2))
            # commands_dict['function'] = self.show_connected_widgets
            return command_list
        elif function == 'set_status':
            arg_list = commands_dict['args']
            command_list = []
            for i, args in enumerate(arg_list):
                command_list.append(lambda arg1=i, arg2=args:
                                    self.set_status(arg1, arg2))
            # commands_dict['function'] = self.show_connected_widgets
            return command_list
        else:
            # print(function)
            raise NotImplementedError

    def widget_connector(self, widget_id, arg_list, func1, func2=None,
                         kwargs1=None, kwargs2=None):
        # var = self.vars[widget_id].get()
        if hasattr(widget_id, 'widget'):
            widget = widget_id.widget
            for i, item in enumerate(self.widgets):
                if widget is item:
                    widget_id = i
                    break
        selected_id = self.widgets[widget_id].current()
        grid_list = arg_list[selected_id]
        show_list = grid_list[0]
        hide_list = grid_list[1]
        self.call_widgets_methods(show_list, func1, kwargs1)
        self.call_widgets_methods(hide_list, func2, kwargs2)

    def show_connected_widgets(self, widget_id, grid_list):
        self.widget_connector(widget_id, grid_list, 'grid', 'grid_remove')

    def set_status(self, widget_id, grid_list):
        self.widget_connector(widget_id, grid_list, 'config', 'config',
                              kwargs1={'state': 'normal'},
                              kwargs2={'state': 'disable'})

    def set_values(self, values, index=0):
        if isinstance(values, (list, tuple, np.ndarray)):
            for i in range(len(self.widgets)):
                if i > len(values) - 1 or i > len(self.widgets) - 1:
                    break
                if values[i] in self.values[i]:
                    self.widgets[i].current(self.values[i].index(values[i]))
                else:
                    raise ValueError('value not found in Combobox')
        else:
            if values in self.values[index]:
                self.widgets[index].current(self.values[index].index(values))
            else:
                raise ValueError('value not found in Combobox')


class EntryButtonSet(MultiEntrySet):
    def __init__(self, frame, label, button_dict, number=1,
                 value=None, **kwargs):
        super().__init__(frame, label, number=number, value=value,
                         justify='left', **kwargs)
        self.dtype = kwargs.pop('dtype', 'string')
        button_factory = button.ButtonFactory()
        self.button = button_factory.create(frame, entry=self.widgets[0],
                                            **button_dict)
        # self.dimensions.grid(row=self.row, column=self.column + number + 1,
        #                      padx=kwargs.get('padx', self.PADX),
        #                      pady=kwargs.get('pady', self.PADY))
        if len(self.columnspan) < (self.number + 2):
            self.columnspan.append(1)

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        row, column = super().set_grid(row=row, column=column, **kwargs)
        column += self.columnspan[self.number]
        self._set_grid(self.button.button, row=row, column=column,
                       columnspan=self.columnspan[self.number + 1],
                       **kwargs)
        return row, column

    # def set_values(self, value, index=0):
    #     if isinstance8
