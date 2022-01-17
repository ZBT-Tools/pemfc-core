import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
import copy

from pemfc.src import global_functions as gf

ID_LIST = []
CONTAINER_LIST = []


def make_list(lst):
    if not isinstance(lst, list):
        out = [lst]
    else:
        out = lst
    return out


def tab_container(child, label, ids):
    tabs = dcc.Tabs(
        [dcc.Tab(child[_], label=label_id, value=tab_id, className='custom-tab',
                 selected_className='custom-tab-selected')
         for _, (label_id, tab_id) in enumerate(zip(label, ids))],
        id='tabs', parent_className='some_container', value=ids[0],
        className='custom-tabs')
    return tabs


def val_container(ids, types='output'):
    inputs = [html.Div(
        [html.Div(id={"type": 'global_children', "id": cont_id},
                  className='gd-desc centered'),
         html.Div(id={"type": 'global_value', "id": cont_id},
                  className='gd-val'),
         html.Div(id={"type": 'global_unit', "id": cont_id},
                  className='gd-un')],
        id={"type": 'global_container', "id": cont_id},
        className='val_container') for cont_id in ids]
    return inputs


def flatten(t):
    return [item for sublist in t for item in sublist]


def spacing(label, dimensions):
    """
    label: xs=33%, s=38%, m=45%, l=50%, xl=70%
    dimensions: s=5%, m=12%
    {n:name, p:percentage}
    """
    spacing_label = {'xs': {'n': 'label-param-s', 'p': 33},
                     's': {'n': 'label-param', 'p': 38},
                     'm': {'n': 'label-param-m', 'p': 45},
                     'l': {'n': 'label-param-l', 'p': 50},
                     'xl': {'n': 'label-param-xl', 'p': 70}, }
    spacing_unit = {'null': {'n': 'ufontnull', 'p': 0},
                    's': {'n': 'ufontm centered', 'p': 5},
                    'm': {'n': 'ufont centered', 'p': 9},
                    'l': {'n': 'ufontl centered', 'p': 12},
                    'xl': {'n': 'ufontxl centered', 'p': 15}}

    p_input = 100 - (spacing_label[label]['p'] + spacing_unit[dimensions]['p'])
    s_input = {'width': f'{p_input}%'}

    return spacing_label[label], spacing_unit[dimensions], s_input


def checklist(checked):
    return [] if checked is False else [1]


# Processing id/val from dict_input to Dash
def id_val_gui_to_dash(label, ids, vals, number, type):
    # id_list = [ids] if not isinstance(ids, list) else ids
    if ids:
        id_list = ['-'.join(ids)] if not isinstance(ids[0], list) else \
            ['-'.join([item for item in inside if isinstance(item, str)]) for
             inside in ids]
    else:
        id_list = []
    dict_ids = {}
    if type == 'EntrySet':
        val_list = [vals] if not isinstance(vals, list) else vals
        num_id = len(id_list)
        num_val = len(val_list)

        if number is None:
            if num_id == num_val:
                number = num_id  # or num_val
            else:
                number = 1

        if number != num_id and num_id == 1:
            id_list = [i_d + f'_{num}' for i_d in id_list
                       for num in range(num_val)]
            # if isinstance(id_list[0], list):
            #     id_list = \
            #         [new_id for sublist in id_list for new_id in sublist
            #          if isinstance(new_id, str)]
            # else:
            #     id_list = [i_d+f'_{num}' for i_d in id_list
            #                    for num in range(number)]

        if number != num_val:
            if isinstance(val_list[0], list):
                val_list = val_list  # multiinput
            elif num_val == 1:
                if num_val != num_id and num_id != number:
                    new_val_list = [val_list * num_id] * num_id
                    if len(flatten(new_val_list)) == number:
                        val_list = new_val_list
                    else:
                        raise ValueError(f'Number of IDs from {label} does not '
                                         f'match with number of value '
                                         f'in the list')
                else:
                    val_list = val_list * num_id  # input

            else:
                raise ValueError(f'Number of IDs from {label} does not match '
                                 f'with number of value in the list')

        # print(label, ' num ', number, ' ids ', len(id_list), ' val ',
        #       len(val_list))
        if len(id_list) == len(val_list):
            if isinstance(val_list[0], list):
                dict_id = [{id_key + f'_{num}': val_value[num] for num, key
                            in enumerate(val_value)} for id_key, val_value
                           in zip(id_list, val_list)]
                [dict_ids.update(d) for d in dict_id]
            else:
                [dict_ids.update({id_l: v_l}) for id_l, v_l in
                 zip(id_list, val_list)]
        else:
            raise IndexError(f'List of ID from {label} does not match with '
                             f'list of value')
    return dict_ids, id_list


# Processing label from list of widget to Dash
def label_gui_to_dash(widget_dicts_list):
    new_widg_dict_l = copy.deepcopy(widget_dicts_list)
    row_check = []
    check_list = []
    for _, widget in enumerate(new_widg_dict_l):
        check_label = widget.get('type', 'EntrySet')
        pos = _
        if check_label == 'Label' and 'row' in widget:
            row = widget.pop('row')
            col = widget.pop('column')
            if col > 1:
                if row == row_check[-1]:
                    pos = check_list[-1]
                    dict_of_widg = new_widg_dict_l[pos]
                    val = make_list(dict_of_widg['label'])
                    val.insert(col - 1, widget['label'])
                    new_widg_dict_l[pos].update({'label': val})
            row_check.append(row)
        check_list.append(pos)

    # unique list to
    unq_list = list(set(check_list))
    num = len(check_list)
    del_list = [x for x in range(num)]
    for ele in sorted(unq_list, reverse=True):
        del del_list[ele]
    for ele in sorted(del_list, reverse=True):
        del new_widg_dict_l[ele]
    return new_widg_dict_l


def implement_widget(kwargs):
    if 'type' in kwargs:
        return row_input(**kwargs)
    else:
        return sub_frame(kwargs)


def frame(tab_dict):
    if 'sub_frame_dicts' in tab_dict:
        return html.Div([sub_frame(subframe) for subframe in
                         tab_dict['sub_frame_dicts']],
                        className='neat-spacing')


def sub_frame(sub_frame_dict):
    bold = 'bolded' if 'bold' in sub_frame_dict.get('font', '') else ''

    if 'highlightbackground' and 'highlightthickness' in sub_frame_dict:
        thickness = sub_frame_dict['highlightthickness']
        bg = sub_frame_dict['highlightbackground']
        border = {'border': f'light{bg} {thickness}px solid'}
        cname = 'box-style'
    else:
        border = {'border': 'initial'}
        cname = None

    if 'sub_frame_dicts' in sub_frame_dict:
        if sub_frame_dict['title'] and sub_frame_dict['show_title'] is True:
            return html.Div(
                [html.Div(children=sub_frame_dict['title'],
                          className=f'title {bold}')] + \
                [html.Div([sub_frame(subframe) for subframe in
                           sub_frame_dict['sub_frame_dicts']])],
                style=border, className=cname)
        else:
            return html.Div([sub_frame(subframe) for subframe in
                             sub_frame_dict['sub_frame_dicts']],
                            style=border, className=cname)

    elif 'widget_dicts' in sub_frame_dict:
        size = {'size_label': sub_frame_dict['size_label']} if 'size_label' \
                in sub_frame_dict else {}
        size_u = {'size_unit': sub_frame_dict['size_unit']} if 'size_unit' \
                  in sub_frame_dict else {}
        new_widg_list = label_gui_to_dash(sub_frame_dict['widget_dicts'])

        # size_u = get_spacing_dimensions(new_widg_list)
        title = sub_frame_dict['title']
        show_title = sub_frame_dict.get('show_title', False)

        spec = sub_frame_dict.get('specifier') if 'specifier' in \
                                                  sub_frame_dict else ''
        id_container = \
            {'type': 'container', 'id': title, 'specifier': spec}

        for widg in new_widg_list:
            widg.update(size) if 'size_label' not in widg else widg.update({})
            widg.update(size_u) if 'size_unit' not in widg else widg.update({})
        # print(new_widg_list)

        if show_title:
            if spec:
                CONTAINER_LIST.append(id_container)
                return html.Div(
                    [html.Div(children=sub_frame_dict['title'],
                              className=f'title {bold}')] + \
                    [implement_widget(widget) for widget in new_widg_list],
                    id=id_container, style=border, className=cname)
                #  implement_widget({**widget, **size, **size_u})
            else:
                return html.Div(
                    [html.Div(children=sub_frame_dict['title'],
                              className=f'title {bold}')] + \
                    [implement_widget(widget) for widget in new_widg_list],
                    style=border, className=cname)
        else:
            if spec:
                CONTAINER_LIST.append(id_container)
                return html.Div(
                    [implement_widget(widget) for widget in new_widg_list],
                    id=id_container, style=border, className=cname)
            return html.Div(
                [implement_widget(widget) for widget in new_widg_list],
                style=border, className=cname)


def get_spacing_dimensions(widget_list):
    space_dimensions = space_label = 0

    for widg in widget_list:
        space_d = widg.get('dimensions', '')
        space_l = widg.get('label', '')
        if space_d:
            space_d_new = len(list(space_d))
            if space_d_new > space_dimensions:
                space_dimensions = space_d_new
        # if space_l:
    if space_dimensions == 0:
        size_unit = {'size_unit': 'null'}
    elif 0 < space_dimensions <= 2:
        size_unit = {'size_unit': 's'}
    elif 2 < space_dimensions <= 6:
        size_unit = {'size_unit': 'm'}
    else:
        size_unit = {'size_unit': 'l'}
    # print(space_dimensions, size_unit)
    return size_unit


def row_input(label='', ids='', value='', type='', dimensions='', options='',
              activated=False, size_label='s', size_unit='s', types='input',
              specifier=False, disabled=False, number=None, **kwargs):
    """
    label: str/list str; ids: str/list of strs for each input ids;
    value: int/float;
    type: checklist, dropdown, or label, or input;
    dimensions: str type, set if unit is used; options:only for dropdown;
    activated: only for checklist, True if checklist is prechecked;
    size_label, size_unit: str type for spacing (see def spacing() );
    types: str type (input/output),id input if not set (trial);
    specifier: id callback more specific, usually widget is set as specifier,
    but more specific specifier can be initialised
    """
    n_ids = kwargs['sim_name'] if 'sim_name' in kwargs else ids
    # print(n_ids)

    dict_ids, id_list = \
        id_val_gui_to_dash(label, n_ids, value, number, type)
    s_label, s_unit, s_input = spacing(size_label, size_unit)

    bold = 'bolded' if 'bold' in kwargs.pop('font', '') else ''

    if specifier == 'visibility':
        visibility = 'visibility'
        specifier = False
    else:
        visibility = ''

    if type == 'EntrySet':  # input
        children = \
            [dbc.Input(
                id={'type': types, 'id': input_id, 'specifier':
                    specifier},
                persistence=True, persistence_type="memory", value=val,
                debounce=True, className='val_input', disabled=disabled)
                for (input_id, val) in
                zip(list(dict_ids.keys()), list(dict_ids.values()))]
    elif type == 'ComboboxSet':  # dropdown
        dd_options = [{'label': val, 'value': val} for val in value] \
            if not options else options
        dd_value = value[0] if not options else value
        children = \
            [dbc.Col(dcc.Dropdown(
                id={'type': types, 'id': input_id, 'specifier': specifier},
                options=dd_options, value=dd_value, persistence=True,
                persistence_type="memory", clearable=False, disabled=disabled,
                className='input-style-dropdown')) for input_id in id_list]
    elif type == 'CheckButtonSet':  # checklist
        value = [value] if not isinstance(value, list) else value
        # print('check val ', checklist(value))
        children = [dbc.Col(dbc.Checklist(
            options=[{"label": "", "value": 1}], value=checklist(val),
            id={"type": types, "id": input_id, 'specifier': specifier},
            inline=True, persistence=True, persistence_type="memory",
            className='checklist')) for input_id, val in zip(id_list, value)]
    elif type == 'Label':  # label
        new_label = [label] if isinstance(label, str) else label
        children = [dbc.Col(html.Div(lbl), className='sm-label') for lbl in
                    new_label]
    else:
        # print(label, type)
        children = ''
        # raise NotImplementedError('Type of Component not implemented')

    # Input
    if children:
        if types == 'multiinput':
            ID_LIST.extend(
                [{'type': types, 'id': input_id, 'specifier': specifier}
                 for input_id in dict_ids.keys()])
        else:
            ID_LIST.extend(
                [{'type': types, 'id': input_id, 'specifier': specifier}
                 for input_id in id_list])
        if visibility == 'visibility':
            # ID container has to make sure that there's only 1 id and number
            id_container = \
                {'type': 'container', 'id': id_list[0], 'specifier':
                    'visibility'}
            CONTAINER_LIST.append(id_container)
            inputs = html.Div(
                [dbc.Label(label, className=s_label['n'] + ' g-0'),
                 html.Div(children, className='centered r_flex g-0',
                          style=s_input),
                 html.Div(dimensions, className=s_unit['n'] + ' g-0')],
                id=id_container,
                className="row-lay r_flex g-0")
        else:
            if type == 'Label':
                if kwargs.get('sticky', '') == 'WNS':
                    inputs = html.Div(label, className=f'section {bold}')
                else:
                    inputs = html.Div(
                        [dbc.Label(className=s_label['n'] + ' g-0'),
                         html.Div(children, className='centered r_flex g-0',
                                  style=s_input),
                         html.Div(dimensions, className=s_unit['n'] + ' g-0')],
                        className="row-lay r_flex g-0")
            else:
                inputs = html.Div(
                    [dbc.Label(label, className=s_label['n'] + ' g-0'),
                     html.Div(children, className='centered r_flex g-0',
                              style=s_input),
                     html.Div(dimensions, className=s_unit['n'] + ' g-0')],
                    className="row-lay r_flex g-0")
    else:
        inputs = html.Div()

    return inputs



