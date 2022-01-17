# global imports
import tkinter as tk
from tkinter import filedialog

# local imports
from . import base
from . import data_transfer


class ButtonFactory:

    @staticmethod
    def create(frame, **kwargs):
        button_type = kwargs.pop('type', None)
        if button_type == 'Button':
            return Button(frame, **kwargs)
        elif button_type == 'RunButton':
            return RunButton(frame, **kwargs)
        elif button_type == 'OpenDirectoryButton':
            return OpenDirectoryButton(frame, kwargs.pop('entry', None),
                                       **kwargs)
        elif button_type == 'OpenFileButton':
            return OpenFileButton(frame, **kwargs)
        elif button_type == 'SaveFileButton':
            return SaveFileButton(frame, **kwargs)
        else:
            raise NotImplementedError('type of WidgetSet not implemented')


class Button(base.Base):

    PADX = 1
    PADY = 1

    def __init__(self, frame, **kwargs):

        label = kwargs.pop('label', '')
        super().__init__(frame, label.lower(), **kwargs)
        self.frame = frame
        kwargs = self.remove_dict_entries(kwargs, self.REMOVE_ARGS)
        self.button = tk.Button(self.frame, text=label, command=self.command,
                                **kwargs)

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        self._set_grid(self.button, row=row, column=column, **kwargs)
        return row, column

    def command(self, *args):
        pass

    def _get_values(self, get_object=False):
        values = {'gui_name': self.name}
        if self.sim_name is not None:
            values['sim_name'] = self.sim_name
        if get_object:
            values['object'] = self
        return values

    def get_values(self, get_object=False):
        return self._get_values(get_object=get_object)


class RunButton(Button):
    def __init__(self, frame, **kwargs):
        super().__init__(frame, **kwargs)


class OpenDirectoryButton(Button):
    def __init__(self, frame, entry=None, **kwargs):
        self.entry = entry
        self.directory = kwargs.pop('directory', None)
        super().__init__(frame, **kwargs)

    def command(self):
        directory = filedialog.askdirectory()
        if directory is None:
            return
        else:
            if isinstance(self.entry, tk.Entry):
                self.entry.delete(0, tk.END)
                self.entry.insert(0, directory)
            return directory


class OpenFileButton(Button):
    def __init__(self, frame, **kwargs):
        self.title = kwargs.pop('title', None)
        self.filetypes = kwargs.pop('filetypes', [])
        super().__init__(frame, **kwargs)

    def command(self):
        try:
            file = filedialog.askopenfile(title=self.title,
                                          filetypes=self.filetypes)
        except AttributeError:
            return
        else:
            return file


class SaveFileButton(Button):
    def __init__(self, frame, file_content=None, **kwargs):
        self.file_content = file_content
        super().__init__(frame, **kwargs)

    def command(self):
        directory = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
        if directory is None:
            return
        else:
            directory.write(self.file_content)
            return directory
