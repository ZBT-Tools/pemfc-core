# import copy
import pathlib
import re
import dash
from dash.dependencies import Input, Output, State, ALL  # ClientsideFunction
from dash import dcc
from dash import html
from dash import dash_table as dt

import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import numpy as np
from pemfc.src import interpolation as ip
from flask_caching import Cache

from pemfc.gui import data_transfer
from pemfc.data import input_dicts
from pemfc import main_app
from pemfc.dash import dash_modal as dm
from pemfc.dash import dash_functions as df
from pemfc.dash import dash_layout as dl

from pemfc.dash.dash_app import app
from dash_tabs import tab1, tab2, tab3, tab4, tab5, tab6

from pemfc import gui

import json
import collections
tabs_list = [tab1.tab_layout, tab2.tab_layout, tab3.tab_layout,
             tab4.tab_layout, tab5.tab_layout, tab6.tab_layout]

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()
# app = dash.Dash(__name__)

# server = app.server


# Setup caching
CACHE_CONFIG = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "simple",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)
# cache2 = Cache()
# cache2.init_app(app.server, config=CACHE_CONFIG)

app.layout = html.Div(
    [html.Div(  # HEADER
        [html.Div(
            html.Div(html.Img(
                    src=app.get_asset_url("logo-zbt-duisburg.png"),
                    id="zbt-image",
                    style={  # "min-height": "60px",
                           "height": "auto",  # "60px",
                           "object-fit": 'contain',
                           'position': 'relative',
                           "width": "auto",
                           "margin": "auto"
                    }),
                id="logoContainer", className="pretty_container"),
            # style={'border': '1px solid grey'},
            className="four columns"),
         html.Div(
             html.Div(html.H3("Fuel Cell Stack Model",
                              style={"margin": "auto",
                                     "min-height": "47px",
                                     # "height": "auto",  # "60px",
                                     "width": "auto",}),
                      className="pretty_container", id="title"),
             # style={'border': '1px solid grey'},
             className="eight columns",)],
        id="header",
        style={  # 'border': '1px solid blue',
               'justify-content': 'space-evenly'}),

     # dcc.Loading(dcc.Store(id="ret_data"), fullscreen=True,
     #             style={"backgroundColor": "transparent"}, type='circle',
     #             color="#0a60c2"),
     dbc.Spinner(dcc.Store(id="ret_data"), fullscreen=True,
                 spinner_style={"width": "10rem", "height": "10rem",
                                "color": "#0a60c2"}),
                 # color="#0a60c2"),
     # empty Div to trigger javascript file for graph resizing
     html.Div(id="output-clientside"),
     # modal for any warning
     dm.modal_axes,

     html.Div(  # MIDDLE
         [html.Div(  # LEFT MIDDLE
             [html.Div(      # LEFT MIDDLE TOP
                 [html.Div(
                     html.Button('Run Simulation', id='run_button',
                                 style={'font-size': '11px', 'margin': 'auto'}),
                     style={'width': '40% !important', 'display': 'flex',
                            'align-items': 'center'}),
                     # dcc.Loading(
                     #     id="run_loading",
                     #     type="circle",
                     #     fullscreen=True,
                     #     children=html.Div(id="run_button")
                     # ),
                  html.Div(
                      [dcc.Dropdown(id='results_dropdown',
                                    className='dropdown_input centered'),
                       dcc.Dropdown(id='results_dropdown_2',
                                    className='dropdown_input centered')],
                      style={'min-width': '60% ', 'min-height': '80px',
                             'display': 'flex', 'flex-direction': 'column',
                             'justify-content': 'space-around'}), ],
                 id='cross-filter-options-new',
                 className='another_pretty_container flex-display'),
              html.Div(  # LEFT MIDDLE MIDDLE
                  [dl.tab_container(
                          tabs_list, label=
                          [k['title'] for k in gui.input.main_frame_dicts],
                          ids=[f'tab{num + 1}' for num in
                               range(len(gui.input.main_frame_dicts))])],

              id='setting_container'),
              html.Div(   # LEFT MIDDLE BOTTOM
                  [html.Div(
                       [html.Div('Settings', className='title'),
                        html.Div(
                           [html.Button('Load', id='load-button'),
                            html.Button('Save', id='save-button')],
                            style={'display': 'flex', 'margin': '5px',
                                   'justify-content': 'space-around'}),
                        dcc.Download(id="savefile-json"),
                        html.Div(dbc.Collapse([
                            dcc.Upload(
                                id='upload-file',
                                children=html.Div([
                                    'Drag and Drop or ',
                                    html.A('Select Files',
                                           style={
                                               'font-weight': 'bold',
                                               'text-decoration': 'underline'})
                                ]), style={'width': '100%',
                                           'height': '60px',
                                           'lineHeight': '60px',
                                           'borderWidth': '1px',
                                           'borderStyle': 'dashed',
                                           'borderRadius': '5px',
                                           'textAlign': 'center',
                                           'margin': 'auto'},
                                # accept='.json',
                                className='dragndrop'),
                            html.Div(id="loading-output-1",
                                     className='output-loading')],
                            id='collapse', is_open=False),)],
                       className='neat-spacing')],
                  className='pretty_container')],
             id="left-column", className='four columns'),

          html.Div(  # RIGHT MIDDLE
              [html.Div(  # RIGHT MIDDLE TOP
                  dl.val_container(
                      ids=['gd1', 'gd2', 'gd3', 'gd4', 'gd5', 'gd6', 'gd7',
                           'gd8', 'gd9', 'gd10']),
                  id="global-data", className='flex-display', ),

               html.Div(  # RIGHT MIDDLE BOTTOM
                   dcc.Graph(id="heatmap_graph"),
                   id='countGraphContainer',
                   className='graph pretty_container'),

               html.Div(
                 [html.Div(
                      [html.Div(
                          [dcc.Dropdown(id='dropdown_line',
                                        placeholder='Choose Plots',
                                        className='dropdown_input'),
                           html.Div(dcc.Dropdown(id='dropdown_line2',
                                                 className='dropdown_input'),
                                    id='dline_div',
                                    style={'visibility': 'hidden'})],
                          style={'margin-bottom': '10px'}),
                       html.Div(
                           dcc.Checklist(id='disp_data',
                                         style={'overflow': 'auto'}),
                           className='display_checklist'),
                       dcc.Store(id='disp_chosen'),
                       dcc.Store(id='disp_clicked'),
                       dcc.Store(id='append_check'),
                       html.Div(
                           [html.Button('Clear List', id='clear_button'),
                            html.Button('Export Data to Table', id='export_b')],
                           style={'display': 'flex', 'flex-direction': 'column',
                                  'margin-bottom': '15px'}),
                       html.Div(
                           [html.Button('Append New Data to Table',
                                        id='append_b'),
                            html.Button('Clear Table', id='clear_table_b')],
                           style={
                               'display': 'flex', 'flex-direction': 'column',
                               'margin-bottom': '5px'})],
                      style={'display': 'flex', 'flex-direction': 'column',
                             'flex': 1}),
                  dcc.Store(id='cells_data'),
                  html.Div(dcc.Graph(id='line_graph'),
                           style={'flex': 4, 'flex-direction': 'column'})],
                 className="pretty_container",
                 style={'display': 'flex',  'flex': '1',
                        'justify-content': 'space-evenly'}), ],
              className="eight columns", id='right-column', )],

         className="flex-display",
         style={'flex-direction': 'row',
                # 'border': '1px solid black',
                'justify-content': 'space-evenly'}, ),
     html.Div(
         [html.Div(dt.DataTable(
             id='table', editable=True, column_selectable='multi'),
             # columns=[{'filter_options': 'sensitive'}]),
             id='div_table', style={'overflow': 'auto'},
             className='pretty_container')],
     style={'position': 'relative', 'margin': '0 0.05% 0 0.7%'})],
    id="mainContainer",
    style={'padding': '5px'})


@cache.memoize()
def simulation_store(**kwargs):
    data_transfer.gui_to_sim_transfer(kwargs, input_dicts.sim_dict)
    global_data, local_data, sim = main_app.main()
    return [global_data[0], local_data[0]]


@app.callback(
    Output('ret_data', 'data'),
    # Output('run_button', 'n_clicks'),
    Input("run_button", "n_clicks"),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'id'),
    State({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'id')

)
def compute_simulation(n_click, inputs, inputs2, ids, ids2):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'run_button' in changed_id and n_click is not None:
        dict_data = df.process_inputs(inputs, inputs2, ids, ids2)
        # print(dict_data)
        datas = {}
        for k, v in dict_data.items():
            datas[k] = {'sim_name': k.split('-'), 'value': v}
        simulation_store(**datas)
    else:
        raise PreventUpdate
    return datas


@app.callback(
    Output({'type': 'input', 'id': ALL, 'specifier': ALL}, 'value'),
    Output({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'value'),
    Output('upload-file', 'contents'),
    Input('upload-file', 'contents'),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'id'),
    State({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'id')
)
def upload_settings(contents, state1, state2, ids, ids2):
    if contents is None:
        raise PreventUpdate
    else:
        try:
            j_file = df.parse_contents(contents)

            list_ids = [id_l['id'] for id_l in ids]
            list_ids2 = [id_l['id'] for id_l in ids2]

            dict_ids = {id_l: num for num, id_l in enumerate(list_ids)}
            dict_ids2 = {id_l: num for num, id_l in enumerate(list_ids2)}

            id_match = set.union(set(list_ids),
                                 set([item[:-2] for item in list_ids2]))

            for k, v in j_file.items():
                if k in id_match:
                    if isinstance(v, list):
                        for num, val in enumerate(v):
                            dict_ids2[k+f'_{num}'] = val
                    else:
                        if isinstance(v, bool):
                            if v is True:
                                dict_ids[k] = [1]
                            else:
                                dict_ids[k] = []
                        else:
                            dict_ids[k] = v
                else:
                    continue

            return list(dict_ids.values()), list(dict_ids2.values()), None
        except Exception as e:
            print(e)
            return state1, state2, None


@app.callback(
    Output("savefile-json", "data"),
    Output('save-button', "n_clicks"),
    Input('save-button', "n_clicks"),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'multiinput', 'id': ALL, 'specifier': ALL}, 'value'),
    State({'type': 'input', 'id': ALL, 'specifier': ALL}, 'id'),
    State({'type': 'multiinput', 'id': ALL, 'specifier':  ALL}, 'id'),
    prevent_initial_call=True,
)
def save_settings(n_clicks, val1, val2, ids, ids2):
    if n_clicks is not None:
        dict_data = df.process_inputs(val1, val2, ids, ids2)  # values first
        sep_id_list = [joined_id.split('-') for joined_id in
                       dict_data.keys()]

        val_list = dict_data.values()
        new_dict = {}
        for sep_id, vals in zip(sep_id_list, val_list):
            current_level = new_dict
            for id_l in sep_id:
                if id_l not in current_level:
                    if id_l != sep_id[-1]:
                        current_level[id_l] = {}
                    else:
                        current_level[id_l] = vals
                current_level = current_level[id_l]

        return dict(content=json.dumps(new_dict, sort_keys=True, indent=2),
                    filename="settings.json"), None


@app.callback(
    [Output({'type': 'global_children', 'id': ALL}, 'children'),
     Output({'type': 'global_value', 'id': ALL}, 'children'),
     Output({'type': 'global_unit', 'id': ALL}, 'children'),
     Output({'type': 'global_container', 'id': ALL}, 'style')],
    Input('ret_data', 'data')
)
def global_outputs(data):
    results = simulation_store(**data)
    g = results[0]
    glob = list(results[0])
    glob_len = len(glob)

    desc = [gl for gl in glob]
    unit = [g[gl]['units'] for gl in glob]
    val = \
        ['{:g}'.format(float('{:.5g}'.format(g[gl]['value'])))
         for gl in glob]
    disp = {'display': 'initial'}
    disps = [{k: v} for k, v in disp.items() for _ in range(glob_len)]
    return desc, val, unit, disps


@app.callback(
    [Output('results_dropdown', 'options'),
     Output('results_dropdown', 'value'),
     Output('dropdown_line', 'options')],
    Input('ret_data', 'data')
)
def get_dropdown_options(data):
    results = simulation_store(**data)
    local_data = results[1]
    values = [{'label': key, 'value': key} for key in local_data if key not in
              ["Channel Location", "Cells", "Cathode",
               "Coolant Channels", "Normalized Flow Distribution"]]
    return values, 'Current Density', values


@app.callback(
    [Output('results_dropdown_2', 'options'),
     Output('results_dropdown_2', 'value')],
    [Input('results_dropdown', 'value'),
     Input('ret_data', 'data')]
)
def get_dropdown_options_2(dropdown_key, data):
    if dropdown_key is None:
        raise PreventUpdate
    else:
        results = simulation_store(**data)
        local_data = results[1]
        if 'value' in local_data[dropdown_key]:
            return [], None
        else:
            options = [{'label': key, 'value': key} for key in
                       local_data[dropdown_key]]
            value = options[0]['value']
            return options, value


@app.callback(
    [Output('dropdown_line2', 'options'),
     Output('dropdown_line2', 'value'),
     Output('dline_div', 'style')],
    [Input('dropdown_line', 'value'),
     Input('ret_data', 'data')]
)
def dropdown_line2(dropdown_key, data):
    if dropdown_key is None:
        raise PreventUpdate
    else:
        results = simulation_store(**data)
        local_data = results[1]
        if 'value' in local_data[dropdown_key]:
            return [], None,  {'visibility': 'hidden'}
        else:
            options = [{'label': key, 'value': key} for key in
                       local_data[dropdown_key]]
            value = options[0]['value']
            return options, value, {'visibility': 'visible'}


@app.callback(
    Output("collapse", "is_open"),
    Input("load-button", "n_clicks"),
    State("collapse", "is_open"),
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    [Output('line_graph', 'figure'),
     Output('cells_data', 'data'),
     Output('disp_data', 'options'),
     Output('disp_data', 'value'),
     Output('disp_chosen', 'data')],
    [Input('dropdown_line', 'value'),
     Input('dropdown_line2', 'value'),
     Input('disp_data', 'value'),
     Input('clear_button', 'n_clicks'),
     Input('line_graph', 'restyleData')],
    [State('ret_data', 'data'),
     State('cells_data', 'data'),
     State('disp_data', 'value'),
     State('disp_chosen', 'data')]
)
def update_line_graph(drop1, drop2, checklist, n_click, rdata,
                      data, state2, state3, state4):
    ctx = dash.callback_context.triggered[0]['prop_id']
    if drop1 is None:
        raise PreventUpdate
    else:

        results = simulation_store(**data)
        local_data = results[1]

        x_key = 'Channel Location'
        xvalues = ip.interpolate_1d(local_data[x_key]['value'])

        if drop1 is None:
            yvalues = np.zeros(len(xvalues))
        else:
            if 'value' in local_data[drop1]:
                yvalues = local_data[drop1]['value']
            elif drop2 is not None:
                yvalues = local_data[drop1][drop2]['value']
            else:
                yvalues = np.zeros((len(local_data['Cell Voltage']['value']),
                                    len(xvalues)))

        fig = go.Figure()
        cells = {}
        for num, yval in enumerate(yvalues):
            fig.add_trace(go.Scatter(x=xvalues, y=yval,
                                     mode='lines+markers',
                                     name=f'Cell {num+1}'))
            cells.update({num+1: {'name': f'Cell {num+1}', 'data': yval}})

        fig.update_layout(xaxis_title='Channel Location [m]',
                          yaxis_title='{} - {}'.format(drop1, drop2))

        options = [{'label': cells[k]['name'], 'value': cells[k]['name']}
                   for k in cells]
        val = sorted([k for k in cells])
        value = [f'Cell {num}' for num in val]
        check = [] if state4 is None else state4

        if 'clear_button.n_clicks' in ctx:
            fig.data = []
            return fig, {}, [], [], []
        else:
            if state3 is None:
                return fig, cells, options, value, check
            else:
                if 'disp_data.value' in ctx:
                    fig.for_each_trace(
                        lambda trace: trace.update(
                            visible=True) if trace.name in state3
                        else trace.update(visible='legendonly'))
                    new_check = list(k['value'] for k in options)
                    [new_check.remove(cell) for cell in state3 if cell in
                     new_check]
                    return fig, cells, options, state3, new_check
                elif 'line_graph.restyleData' in ctx:
                    read = rdata[0]['visible']
                    read_num = rdata[1][0]

                    if len(read) == 1:
                        if isinstance(read[0], str):  # lose (legendonly)
                            if f'Cell {read_num+1}' not in check:
                                check.append(f'Cell {read_num+1}')
                        else:  # isinstance(read, bool): #add (True)
                            try:
                                if 'Cell {}'.format(read_num + 1) in check:
                                    check.remove('Cell {}'.format(read_num + 1))
                            except ValueError:
                                pass
                        [value.remove(val) for val in check if val in value]
                        fig.for_each_trace(
                            lambda trace: trace.update(
                                visible=True) if trace.name in value
                            else trace.update(visible='legendonly'))

                        return fig, cells, options, value, check
                    else:
                        check_new = [f'Cell {x[0]+1}' for x in enumerate(read)
                                     if x[1] == 'legendonly']
                        [value.remove(che) for che in check_new
                         if che in value]
                        fig.for_each_trace(
                            lambda trace: trace.update(
                                visible=True) if trace.name in value
                            else trace.update(visible='legendonly'))

                        return fig, cells, options, value, check_new
                else:
                    return fig, cells, options, value, check


@app.callback(
    [Output('table', 'columns'),
     Output('table', 'data'),
     Output('table', 'export_format'),
     Output('append_check', 'data')],
    [Input('disp_data', 'value'),  # from display
     Input('cells_data', 'data'),  # from line graph
     Input('disp_chosen', 'data'),  # from line graph (but only the val of cell)
     Input('export_b', 'n_clicks'),  # button1
     Input('append_b', 'n_clicks'),  # button2
     Input('clear_table_b', 'n_clicks')],
    [State('ret_data', 'data'),
     State('table', 'columns'),
     State('table', 'data'),
     State('append_check', 'data')]
)
def list_to_table(val, data, data2, n1, n2, n3, state, state2, state3, state4):
    ctx = dash.callback_context.triggered[0]['prop_id']
    if val is None:
        raise PreventUpdate
    else:
        results = simulation_store(**state)
        local_data = results[1]
        digit_list = \
            sorted([int(re.sub('[^0-9\.]', '', inside)) for inside in val])

        index = [{'id': 'Channel Location', 'name': 'Channel Location',
                  'deletable': True}]
        columns = [{'deletable': True, 'renamable': True,
                    'selectable': True, 'name': 'Cell {}'.format(d),
                    'id': 'Cell {}'.format(d)} for d in digit_list]
        # list with nested dict

        x_key = 'Channel Location'
        xvalues = ip.interpolate_1d(local_data[x_key]['value'])
        datas = [{**{'Channel Location': cell},
                  **{data[k]['name']: data[k]['data'][num] for k in data}}
                 for num, cell in enumerate(xvalues)]  # list with nested dict

        if state4 is None:
            appended = 0
        else:
            appended = state4

        if 'export_b.n_clicks' in ctx:
            return index+columns, datas, 'csv', appended
        elif 'clear_table_b.n_clicks' in ctx:
            return [], [], 'none', appended
        elif 'append_b.n_clicks' in ctx:
            if n1 is None or state3 == [] or state2 == [] or \
                    ctx == 'clear_table_b.n_clicks':
                raise PreventUpdate
            else:
                appended += 1
                app_columns = \
                    [{'deletable': True, 'renamable': True,
                      'selectable': True, 'name': 'Cell {}'.format(d),
                      'id': 'Cell {}'.format(d) + '-'+str(appended)}
                     for d in digit_list]
                new_columns = state2 + app_columns
                app_datas = \
                    [{**{'Channel Location': cell},
                      **{data[k]['name']+'-'+str(appended): data[k]['data'][num]
                         for k in data}}
                     for num, cell in enumerate(xvalues)]
                new_data_list = []
                new_datas = \
                    [{**state3[i], **app_datas[i]}
                     if state3[i][x_key] == app_datas[i][x_key]
                     else new_data_list.extend([state3[i], app_datas[i]])
                     for i in range(len(app_datas))]
                new_datas = list(filter(None.__ne__, new_datas))
                new_datas.extend(new_data_list)

                return new_columns, new_datas, 'csv', appended
        else:
            if n1 is None or state2 == []:
                return state2, state3, 'none', appended
            else:
                return state2, state3, 'csv', appended


@app.callback(
    Output("heatmap_graph", "figure"),
    [Input('results_dropdown', 'value'), Input('results_dropdown_2', 'value'),
     Input('ret_data', 'data')],
)
def update_graph(dropdown_key, dropdown_key_2, data):
    if dropdown_key is None:
        raise PreventUpdate
    else:
        results = simulation_store(**data)
        local_data = results[1]
        x_key = 'Channel Location'
        y_key = 'Cells'
        xvalues = ip.interpolate_1d(local_data[x_key]['value'])
        yvalues = local_data[y_key]['value']
        if dropdown_key is None:
            zvalues = np.zeros((len(xvalues), len(yvalues)))
        else:
            if 'value' in local_data[dropdown_key]:
                zvalues = local_data[dropdown_key]['value']
            elif dropdown_key_2 is not None:
                zvalues = local_data[dropdown_key][dropdown_key_2]['value']
            else:
                zvalues = np.zeros((len(xvalues), len(yvalues)))
            # else:
            #     zvalues = local_data[dropdown_key][dropdown_key_2]['value']

        fig = go.Figure(go.Heatmap(z=zvalues, x=xvalues, y=yvalues,
                                   xgap=1, ygap=1))
        fig.update_xaxes(showgrid=True, tickmode='array',
                         tickvals=local_data[x_key]['value'])
        fig.update_yaxes(showgrid=True, tickmode='array', tickvals=yvalues)
    return fig


if __name__ == "__main__":
    # [print(num, x) for num, x in enumerate(dl.ID_LIST) ]
    app.run_server(debug=True, use_reloader=False)
