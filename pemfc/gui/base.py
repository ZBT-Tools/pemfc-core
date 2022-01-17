# global imports
from abc import ABC
from tkinter import Grid
import copy

# local imports
from pemfc import global_functions as gf


class Base(ABC):

    PADX = 1
    PADY = 1

    REMOVE_ARGS = ['row', 'column', 'grid_location', 'columnspan',
                   'rowspan', 'sticky', 'sim_name', 'dtype', 'padx', 'pady',
                   'width', 'weights', 'size_label', 'size_unit', 'specifier',
                   'types']

    def __init__(self, master, name, **kwargs):
        if hasattr(master, 'name'):
            self.name = copy.deepcopy(gf.ensure_list(master.name))
            self.name.append(name.lower())
        else:
            self.name = [name.lower()]
        self.sim_name = kwargs.pop('sim_name', None)
        self.padx = kwargs.pop('padx', self.PADX)
        self.pady = kwargs.pop('pady', self.PADY)
        self.rowspan = kwargs.pop('rowspan', 1)
        self.columnspan = kwargs.pop('columnspan', 1)
        grid_location = kwargs.pop('grid_location', (None, 0))
        self.row = kwargs.pop('row', grid_location[0])
        self.column = kwargs.pop('column', grid_location[1])
        self.sticky = kwargs.pop('sticky', 'NE')
        self.weights = kwargs.pop('weights', None)

    def _set_grid(self, widget, **kwargs):
        # Grid.rowconfigure(self.frame, self.row, weight=1)
        # Grid.columnconfigure(self.frame, self.column, weight=1)

        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        # self.frame.rowconfigure(row, weight=1)
        # self.frame.columnconfigure(column, weight=1)
        widget.grid(row=row, column=column,
                    padx=kwargs.pop('padx', self.padx),
                    pady=kwargs.pop('pady', self.pady),
                    columnspan=kwargs.pop('columnspan', self.columnspan),
                    rowspan=kwargs.pop('rowspan', self.rowspan),
                    sticky=kwargs.pop('sticky', self.sticky), **kwargs)
        return row, column

    @staticmethod
    def remove_dict_entries(dictionary, entries):
        for entry in entries:
            dictionary.pop(entry, None)
        return dictionary

    def call_commands(self):
        pass

