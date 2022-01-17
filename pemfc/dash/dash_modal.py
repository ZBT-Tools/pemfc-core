import dash
from dash.dependencies import Input, Output, State
from dash import html
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

from pemfc.dash.dash_app import app


modal_axes = html.Div([
    # Progress Modal
    dbc.Modal(
        [dbc.ModalHeader(dbc.ModalTitle("Work in Progress!",
                                        style={'font-weight': 'bold'})),
         dbc.ModalBody([
             html.Div("Data from JSON file has been loaded!"),
             html.Div("You can click 'Run Simulation' to simulate the loaded "
                      "parameter or continue to change the loaded parameter "
                      "value")
         ]),
         dbc.ModalFooter(
             dbc.Button("Close", id="progress-close",
                        className="ms-auto", n_clicks=0)), ],
        id="progress-modal", is_open=False),

    # Load Error Modal
    dbc.Modal(
        [dbc.ModalHeader(dbc.ModalTitle("Attention!", style={'font-weight':
                                                                 'bold'})),
         dbc.ModalBody("Upload function only accepts JSON file!"),
         dbc.ModalFooter(
             dbc.Button("Close", id="load-error-close",
                        className="ms-auto", n_clicks=0)), ],
        id="load-error-modal", is_open=False),
    # Load wrong File
    dbc.Modal(
        [dbc.ModalHeader(dbc.ModalTitle("Attention!",
                                        style={'font-weight': 'bold'})),
         dbc.ModalBody([
             html.Div("IDs from JSON file is not compatible/does not match up"),
             html.Div("Please review the uploaded JSON file again")
         ]),
         dbc.ModalFooter(
             dbc.Button("Close", id="load-wrong-error-close",
                        className="ms-auto", n_clicks=0)), ],
        id="load-wrong-modal", is_open=False),

    # Append Error Modal
    dbc.Modal(
        [dbc.ModalHeader(dbc.ModalTitle("Attention!", style={'font-weight':
                                                             'bold'})),
         dbc.ModalBody([
             html.Div("Data cannot be appended to an empty table."),
             html.Div("Please export the data first to the table!")
         ]),
         dbc.ModalFooter(
             dbc.Button("Close", id="append-error-close",
                        className="ms-auto", n_clicks=0)),],
        id="append-error-modal", is_open=False)])


@app.callback(
    [Output("load-error-modal", "is_open"),
     Output("progress-modal", "is_open")],
    [Input('upload-file', 'contents'),
     Input("load-error-close", "n_clicks"),
     Input("progress-close", "n_clicks")],
    [State('upload-file', 'filename'),
     State("load-error-modal", "is_open"),
     State("progress-modal", "is_open")],
)
def toggle_modal(contents, n1, n2, state, is_open, is_open2):
    if contents is not None:
        if 'json' not in state or n1:
            return not is_open, is_open2
        if 'json' in state or n2:
            return is_open, not is_open2
        return is_open, is_open2
    else:
        raise PreventUpdate


@app.callback(
    Output("append-error-modal", "is_open"),
    [Input('append_b', 'n_clicks'), Input("append-error-close", "n_clicks")],
    [State('table', 'data'), State("append-error-modal", "is_open")],
)
def toggle_modal(appended, n1, state, is_open):
    ctx = dash.callback_context.triggered[0]['prop_id']
    if state is None or state == []:
        if 'append_b.n_clicks' in ctx or n1:
            return not is_open
    return is_open


