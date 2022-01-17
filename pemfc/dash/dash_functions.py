from functools import wraps
from flask_caching import Cache
import base64
import io
import json
import copy
import collections
from glom import glom

from pemfc.data import input_dicts
from pemfc import gui


def unstringify(val):
    """
    Is used to change any str value created by DBC.Input once initialised due to
    not defining the component as type Number.
    """
    if isinstance(val, str):
        if '.' in val:
            try:
                yield float(val)
            except (ValueError, NameError):
                yield val
        else:
            try:
                yield int(val)
            except (ValueError, NameError):
                yield val
    elif isinstance(val, list):
        try:
            yield list((float(v) for v in val))
        except ValueError:
            yield list((v for v in val))
    else:
        yield val


# def multi_inputs(inputs):
#     input_list = []
#     for num, val in enumerate(inputs):
#         if (num % 2) != 0:
#             input_list.append(inputs[num - 1:num + 1])
#     return input_list


def multi_inputs(dicts):
    dict_list = {}
    for k, v in dicts.items():
        if k[-1:].isnumeric() is False:
            dict_list.update({k: v})
        else:
            if k[:-2] not in dict_list:
                dict_list.update({k[:-2]: v})
            elif k[:-2] in dict_list:
                if not isinstance(dict_list[k[:-2]], list):
                    new_list = [dict_list[k[:-2]]]
                else:
                    new_list = dict_list[k[:-2]]
                new_list.append(v)
                dict_list.update({k[:-2]: new_list})
    return dict_list


def dict_inputs(value='', ids=''):
    data = {}
    if value != '':
        for id_l, val in zip(ids, value):
            checked_val = list(unstringify(val))
            data[id_l['id']] = checked_val[0]
    else:
        for val, id_l in enumerate(ids):
            data[id_l['id']] = val
    return data


def list_dict_inputs(value='', ids=''):
    data = []
    if value != '':
        for id_l, val in zip(ids, value):
            checked_val = list(unstringify(val))
            data.append({id_l['id']: checked_val[0]})
    else:
        for val, id_l in enumerate(ids):
            data.append({id_l['id']: val})
    return data


def multival_input(value, input_ids):
    inputs = [list(unstringify(v))[0] for v in value]
    data_dict = {ids['id']: inp for ids, inp in zip(input_ids, inputs)}
    set_list = list(set([k[:-2] for k, v in data_dict.items()]))
    data = {k: [data_dict[k + '-z'], data_dict[k + '-x']] for k in set_list}

    return {k: {'sim_name': k.split('-'), 'value': v} for k, v in data.items()}


def dash_kwarg(inputs):
    """
    :param inputs:
    :return:
    """
    def accept_func(func):
        @wraps(func)
        def wrapper(*args):
            input_names = [item.component_id for item in inputs]
            kwargs_dict = dict(zip(input_names, args))
            return func(**kwargs_dict)
        return wrapper
    return accept_func


# def tab_compile(inputs, ):


def compile_data(**kwargs):
    """
    Used in compiling data from each components into a dcc.Store for each Tab
    """
    dash_dict = collections.defaultdict(dict)
    for k, v in kwargs.items():
        t = dash_dict[k]['sim_name'] = k.split("-")
        c_v = list(unstringify(v))
        if t[-1] in ['cool_flow', 'cool_ch_bc', 'calc_distribution']:
            dash_dict[k]['value'] = bool(c_v[0])
        else:
            dash_dict[k]['value'] = c_v[0]
    return dict(dash_dict)


# def extract_component(widget, source):
#     if widget not in source:


def gen_dict_extract(key, var):
    if hasattr(var, 'items'):
        for k, v in var.items():
            if k == key:
                yield var
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result
    elif isinstance(var, (list, tuple)):
        for d in var:
            for result in gen_dict_extract(key, d):
                yield result


def parse_contents(contents):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    j_file = json.load(io.StringIO(decoded.decode('utf-8')))
    # val_j_file =
    # [print(i, x) for i, x in enumerate(j_file)]

    source_dict = copy.deepcopy(gui.input.main_frame_dicts)
    target_dict = copy.deepcopy(input_dicts.sim_dict)

    settings, name_lists = \
        gui.data_transfer.gui_to_sim_transfer(source_dict, target_dict)

    js_out = {}
    for n in name_lists:
        try:
            js_out.update({'-'.join(n): glom(j_file, '.'.join(n))})
        except Exception:
            pass

    return js_out


def process_inputs(inputs, multiinputs, id_inputs, id_multiinputs):
    new_inputs = []
    for val in inputs + multiinputs:
        new_val = list(unstringify(val))[0]

        if isinstance(new_val, list):
            if len(new_val) == 0:
                new_val = bool(new_val)
            else:
                if len(new_val) == 1 and new_val[0] == 1:
                    new_val = bool(new_val)
        new_inputs.append(new_val)
    # print(new_inputs)

    new_ids = [id_l['id'] for id_l in id_inputs] + \
              [id_l['id'] for id_l in id_multiinputs]
    # [print(x) for x in new_ids]
    dict_data = {}
    for id_l, v_l in zip(new_ids, new_inputs):
        dict_data.update({id_l: v_l})
    new_dict_data = multi_inputs(dict_data)
    return new_dict_data



# def parse_contents(contents):
#     content_type, content_string = contents.split(',')
#
#     decoded = base64.b64decode(content_string)
#     try:
#         j_file = json.load(io.StringIO(decoded.decode('utf-8')))
#         source_dict = copy.deepcopy(gui.input.main_frame_dicts)
#         target_dict = copy.deepcopy(input_dicts.sim_dict)
#
#         settings, name_lists = \
#             gui.data_transfer.gui_to_sim_transfer(source_dict, target_dict)
#         js_out = {'-'.join(n): glom(j_file, '.'.join(n)) for
#                   n in name_lists if 'loss' not in n[-1] and 'directory' not
#                   in n[-1] and 'calc_current_density' not in n[-1] and
#                   'calc_temperature' not in n[-1]}
#         return js_out
#     except Exception as e:
#         return f'Error {e}'
        # return f'{e}: There was an error processing this file.'
