# global imports
import tkinter as tk
from tkinter import ttk
import os
import sys
import json

# local imports
from pemfc.gui import frame
from pemfc.gui import input
from pemfc.gui import data_transfer
# from pemfc.gui import icon
from pemfc import main_app
from pemfc.data import input_dicts


class NotebookApp:
    def __init__(self, master, main_frame_dicts=None, **kwargs):
        self.notebook = ttk.Notebook(master)
        self.master = master
        self.registry = {}
        # self.configure_gui()
        # Make tabs and corresponding frames
        self.frame_factory = frame.FrameFactory()
        self.frames = []
        if main_frame_dicts is not None:
            for i, main_frame_dict in enumerate(main_frame_dicts):
                main_frame = self.frame_factory.create_frame(self.notebook,
                                                             **main_frame_dict)
                self.frames.append(main_frame)
                self.notebook.add(main_frame, text=main_frame_dict['title'],
                                  sticky='WENS')

        sim_frame = self.frames[-1]

        # configure load settings button
        self.load_settings_button = \
            sim_frame.sub_frames[0].sub_frames[-2].widgets[0]
        self.load_settings_button.button.configure(command=self.load_settings)
        self.load_settings_button.filetypes = [('JSON Files', ['.json'])]
        self.load_settings_button.title = 'Please select settings file.'

        # configure save settings button
        self.save_settings_button = \
            sim_frame.sub_frames[0].sub_frames[-2].widgets[1]
        self.save_settings_button.button.configure(command=self.save_settings)

        # configure run button
        run_frame = sim_frame.sub_frames[0].sub_frames[-1]
        self.run_button = run_frame.widgets[0]
        self.run_button.button.configure(command=self.run)

        # # add progress bar
        # self.progress = ttk.Progressbar(run_frame, orient=tk.HORIZONTAL,
        #                            length=100, mode='determinate')
        # run_frame.widgets.append(self.progress)

        self.notebook.select(self.frames[0])
        self.notebook.enable_traversal()
        self.set_grid()
        self.call_commands()

        # set custom grid
        # values = self.get_values(get_object=True)
        # widget_set = values['physical properties']['physical properties']\
        #     ['porous layers']['gas diffusion layer porosity']['object']
        # widget_set.widgets[0].grid(column=1, columnspan=2)
        # widget_set.widgets[1].grid(column=2, columnspan=2)

    def set_grid(self, grid_list=None, **kwargs):
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.columnconfigure(0, weight=1)
        # self.notebook.grid(sticky='WENS', **kwargs)
        for fr in self.frames:
            fr.set_grid(grid_list=grid_list, **kwargs)
        self.notebook.grid(sticky='WEN', **kwargs)

    def get_values(self, get_object=False):
        # return [fr.get_values() for fr in self.frames]
        return {fr.name[-1]: fr.get_values(get_object=get_object)
                for fr in self.frames}

    def call_commands(self):
        for item in self.frames:
            item.call_commands()

    # def configure_gui(self):
    #     self.master.tk.call('wm', 'iconphoto', self.master._w,
    #                         tk.PhotoImage(data=icon.icon()))
    #     self.master.title("Example")
    #     self.master.minsize(250, 50)

    def get_settings(self):
        values = self.get_values()
        settings, name_lists = \
            data_transfer.gui_to_sim_transfer(values, input_dicts.sim_dict)
        return settings

    def save_settings(self):
        settings_dict = self.get_settings()
        self.save_settings_button.file_content = \
            json.dumps(settings_dict, indent=2)
        self.save_settings_button.command()
        return settings_dict

    def load_settings(self):
        file_directory = self.load_settings_button.command()
        if file_directory is None:
            return
        else:
            settings_dict = json.load(file_directory)
            widgets_registry = self.get_values(get_object=True)
            data_transfer.sim_to_gui_transfer(settings_dict, widgets_registry)
            self.call_commands()
            input_dicts.sim_dict.update(settings_dict)

    def run(self):
        settings = self.get_settings()
        main_app.main()


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS',
                        os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


if __name__ == "__main__":
    root = tk.Tk()
    root.title('PEMFC Model')
    if getattr(sys, 'frozen', False):
        # frozen
        gui_dir = os.path.dirname(sys.executable)
    else:
        # unfrozen
        gui_dir = os.path.dirname(os.path.abspath(__file__))
    # gui_dir = os.path.dirname(os.path.abspath(__file__))
    root.iconbitmap(os.path.join(gui_dir, 'logo-zbt.ico'))

    # root.resizable(False, False)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    # style = ttk.Style()

    # style.theme_create('style', #parent='alt',
    #                         settings={
    #                             'Combobox': {
    #                                 'configure': {
    #                                     #'fieldbackground': 'white',
    #                                     #'selectbackground': 'white',
    #                                     #'selectforeground': 'black'
    #                                     #'highlightbackground': 'Yellow'
    #                                 }
    #                             }
    #                         }
    #                         )
    #
    # style.theme_use('style')

    base_app = NotebookApp(root, main_frame_dicts=input.main_frame_dicts)

    root.mainloop()
