
from dash import html

from pemfc import gui
from pemfc.dash import dash_layout as dl


tab_layout = html.Div(dl.frame(gui.input.main_frame_dicts[3]))

