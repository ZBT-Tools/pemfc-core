# global imports
import tkinter as tk
from tkinter import ttk
import numpy as np

# local imports
from . import base
from . import widget_set as ws
from . import button
from . import widget_factory as wf


class FrameFactory:
    @staticmethod
    def create_frame(master, sub_frame_dicts: list = None,
                     widget_dicts: list = None,  # button_dicts: list = None,
                     **kwargs):
        # if button_dicts is None:
        if sub_frame_dicts is None:
            return BaseFrame(master, kwargs.pop('title'),
                             widget_dicts=widget_dicts, **kwargs)
        else:
            return MainFrame(master, kwargs.pop('title'),
                             sub_frame_dicts=sub_frame_dicts,
                             widget_dicts=widget_dicts, **kwargs)
        # else:
        #     return ButtonMainFrame(master, widget_dicts=widget_set_dicts,
        #                            button_dicts=button_dicts, **kwargs)


class BaseFrame(base.Base, tk.Frame):

    PADX = 2
    PADY = 2

    def __init__(self, master, name, widget_dicts: list = None, **kwargs):
        show_title = kwargs.pop('show_title', False)
        title = name  # kwargs.pop('title', None)
        name = name.lower()
        font = kwargs.pop('font', None)
        command_order = kwargs.pop('command_order', None)
        self.initialize = True
        if isinstance(master, ttk.Notebook):
            self.notebook_tab = True
        else:
            self.notebook_tab = False
        self.title = None
        base.Base.__init__(self, master, name,
                           sticky=kwargs.pop('sticky', 'WENS'), **kwargs)
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)

        tk.Frame.__init__(self, master, name=name, **kwargs)
        if show_title:
            self.title = self.set_title(title, font=font)

        # self.grid(sticky=kwargs.pop('sticky', self.sticky),
        #           row=kwargs.pop('row', self.grid_location[0]),
        #           column=kwargs.pop('column', self.grid_location[1]),
        #           padx=kwargs.get('padx', self.PADX),
        #           pady=kwargs.get('pady', self.PADY))
        self.widget_factory = wf.WidgetFactory()
        self.widgets = []
        if widget_dicts is not None:
            self.widgets = [self.widget_factory.create(self, **w_dict)
                            for w_dict in widget_dicts]
            # self.columns = np.max([widget.columns for widget in self.widgets])
        self.widget_grid = []
        if command_order is None:
            self.command_order = None
        elif isinstance(command_order, (list, tuple)):
            if len(command_order) == len(self.widgets):
                self.command_order = command_order
            else:
                raise ValueError('argument command_order must have the same '
                                 'number of entries as widgets in this frame')
        else:
            raise TypeError('argument command_order must be of type list')

    def set_title(self, text, font=None, **kwargs):
        title = ws.Label(self, label=text, font=font, sticky='WENS', **kwargs)
        return title

    @staticmethod
    def _get_values(tk_objects, get_object=False):
        values = {}
        for item in tk_objects:
            if hasattr(item, 'get_values'):
                values[item.name[-1]] = item.get_values(get_object=get_object)
        return values

    def get_values(self, tk_objects=None, get_object=False):
        if tk_objects is None:
            tk_objects = self.widgets
        return self._get_values(tk_objects, get_object=get_object)

    def set_grid(self, tk_objects=None, grid_list=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        if row is None:
            row = 0
        if column is None:
            column = 0
        self.rowconfigure(row, weight=1)
        self.columnconfigure(column, weight=1)

        if not self.notebook_tab:
            # self.grid(sticky=kwargs.pop('sticky', self.sticky),
            #           row=row, column=column,
            #           padx=kwargs.get('padx', self.PADX),
            #           pady=kwargs.get('pady', self.PADY), **kwargs)
            self._set_grid(self, **kwargs)
        if tk_objects is None:
            tk_objects = self.widgets
        if grid_list is None:
            if tk_objects and self.title is None:
                row = -1
            else:
                row = 0
            column = 0

            for i, tk_object in enumerate(tk_objects):

                if hasattr(tk_object, 'row') and tk_object.row is not None:
                    row = tk_object.row
                else:
                    row += 1
                if hasattr(tk_object, 'column') \
                        and tk_object.column is not None:
                    column = tk_object.column
                else:
                    column = 0
                tk_object.row = row
                tk_object.column = column
                if hasattr(tk_object, 'set_grid'):
                    tk_object.set_grid(**kwargs)
                else:
                    tk_object.grid(**kwargs)
                if hasattr(tk_object, 'weights') \
                        and tk_object.weights is not None:
                    if isinstance(tk_object.weights, (list, tuple)):
                        self.rowconfigure(row, weight=tk_object.weights[0])
                        self.columnconfigure(column,
                                             weight=tk_object.weights[1])
                # else:
                #     self.rowconfigure(row, weight=1)
                    # self.columnconfigure(column, weight=1)

            if self.title is not None:
                self.title.set_grid(row=0, column=0,
                                    columnspan=self.grid_size()[0])
        # for i in range(self.grid_size()[1]):
        #     for j in range(self.grid_size()[0]):
        #         self.rowconfigure(i, weight=1)
        #         self.columnconfigure(j, weight=1)
        self.add_widget_grid()
        return row, column

    def add_widget_grid(self):
        grid_size = self.grid_size()
        self.widget_grid = [[None for column in range(grid_size[0])]
                            for row in range(grid_size[1])]
        widget_list = self.winfo_children()
        for widget in widget_list:
            row = widget.grid_info()['row']
            column = widget.grid_info()['column']
            self.widget_grid[row][column] = widget

    def add_widget(self, widget):
        self.widgets.append(widget)

    def call_commands(self):
        if self.command_order is not None and self.initialize:
            for i in range(len(self.widgets)):
                j = self.command_order[i]
                widget = self.widgets[j]
                if hasattr(widget, 'call_commands'):
                    widget.call_commands()
        else:
            if self.command_order is None:
                for widget in self.widgets:
                    if hasattr(widget, 'call_commands'):
                        widget.call_commands()
        self.initialize = False


class MainFrame(BaseFrame):
    def __init__(self, master, name, sub_frame_dicts: list = None,
                 widget_dicts: list = None, **kwargs):
        super().__init__(master, name, widget_dicts=widget_dicts, **kwargs)
        self.frame_factory = FrameFactory()
        self.sub_frames = []
        if sub_frame_dicts is not None:
            for i, frame_dict in enumerate(sub_frame_dicts):
                sub_frame = self.frame_factory.create_frame(self, **frame_dict)
                self.sub_frames.append(sub_frame)

    def set_grid(self, grid_list=None, **kwargs):
        row, column = super().set_grid(tk_objects=self.sub_frames,
                                       grid_list=grid_list, **kwargs)
        if self.widgets:
            kwargs.pop('row', None)
            kwargs.pop('column', None)
            row, column = super().set_grid(tk_objects=self.widgets,
                                           grid_list=grid_list,
                                           row=row, column=column, **kwargs)
        return row, column

    def get_values(self, tk_objects=None, get_object=False):
        if tk_objects is None:
            tk_objects = self.sub_frames
        values = self._get_values(tk_objects=tk_objects, get_object=get_object)
        if self.widgets:
            values.update(self._get_values(tk_objects=self.widgets,
                                           get_object=get_object))
        return values

    def call_commands(self):
        for frame in self.sub_frames:
            frame.call_commands()


# class ButtonMainFrame(BaseFrame):
#     def __init__(self, master, sub_frame_dicts: list = None,
#                  widget_dicts: list = None, button_dicts: list = None,
#                  **kwargs):
#         super().__init__(master, widget_dicts=widget_dicts, **kwargs)
#
#         button_factory = button.ButtonFactory()
#         self.buttons = []
#         if button_dicts is not None:
#             self.buttons = [button_factory.create(self, **b_dict)
#                             for b_dict in button_dicts]
#
#     def set_grid(self, grid_list=None, **kwargs):
#         row, column = super().set_grid(grid_list=grid_list, **kwargs)
#         if self.buttons:
#             kwargs.pop('row', None)
#             kwargs.pop('column', None)
#             row, column = BaseFrame.set_grid(self, tk_objects=self.buttons,
#                                              grid_list=grid_list,
#                                              row=row, column=column, **kwargs)
#         return row, column
#
#     def get_values(self, tk_objects=None):
#         values = super().get_values()
#         if self.buttons:
#             values.update(self._get_values(tk_objects=self.buttons))
#         return values

