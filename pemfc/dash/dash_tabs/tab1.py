
from dash.dependencies import Input, Output, State, ALL
from dash import html

from pemfc.dash.dash_app import app
from pemfc import gui
from pemfc.dash import dash_layout as dl

tab_layout = html.Div(dl.frame(gui.input.main_frame_dicts[0]))


@app.callback(
    Output({'type': ALL, 'id': ALL, 'specifier': 'disable_basewidth'},
           'disabled'),
    Input({'type': ALL, 'id': ALL, 'specifier': 'dropdown_activate_basewidth'},
          'value'))
def disabled_callback(value):
    for num, val in enumerate(value):
        if val == "trapezoidal":
            value[num] = False
        else:
            value[num] = True
    return value
