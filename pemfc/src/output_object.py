import string
import weakref
from copy import deepcopy


class OutputObject:

    # PRINT_HIERARCHY = 3
    # CLUSTER_NAMES = [['Cell', 'Flow Circuit']]
    _instances = set()

    def __init__(self, name, **kwargs):
        super().__init__()
        assert isinstance(name, str)
        self._name = name
        self.active = True

        self.single_print_data = {}
        self.multi_print_data = {}
        self.print_data = [self.single_print_data, self.multi_print_data]
        self._instances.add(weakref.ref(self))

    def _get_name(self):
        return self._name

    def _set_name(self, name):
        self._name = name

    @property
    def name(self):
        return self._get_name()

    @name.setter
    def name(self, name):
        self._set_name(name)

    def extend_data_names(self, name, prepend=True):
        for i, print_data in enumerate(self.print_data):
            print_data_new_keys = {}
            for key, value in print_data.items():
                if prepend:
                    new_key = name + ' ' + key
                else:
                    new_key = key + ' ' + name
                print_data_new_keys[new_key] = value
                # print_data.pop(key)
            self.print_data[i] = print_data_new_keys
        # print_data_new_keys = {}
        # for key, value in self.print_data_1d.items():
        #     if prepend:
        #         new_key = name + ' ' + key
        #     else:
        #         new_key = key + ' ' + name
        #     print_data_new_keys[new_key] = value
        # self.print_data_1d = print_data_new_keys
        # print_data_new_keys = {}
        # for key, value in self.print_data_2d.items():
        #     if prepend:
        #         new_key = name + ' ' + key
        #     else:
        #         new_key = key + ' ' + name
        #     print_data_new_keys[new_key] = value
        # self.print_data_2d = print_data_new_keys
        # self.print_data = [self.print_data_1d, self.print_data_2d]

    @classmethod
    def getinstances(cls):
        dead = set()
        for ref in cls._instances:
            obj = ref()
            if obj is not None:
                yield obj
            else:
                dead.add(ref)
        cls._instances -= dead

    def copy(self):
        copy = deepcopy(self)
        self._instances.add(weakref.ref(copy))
        return copy

    def add_print_data(self, data_array, name, units='-', plot_axis=-1,
                       sub_names=None, multi_data=False):
        if sub_names is not None or multi_data:
            if sub_names is None:
                sub_names = [str(i+1) for i in range(len(data_array))]
            self.multi_print_data[name] = \
                {sub_names[i]:
                 {'value': data_array[i], 'units': str(units),
                  'plot_axis': plot_axis, 'save': True}
                 for i in range(len(sub_names))}
        else:
            self.single_print_data[name] = \
                {'value': data_array, 'units': str(units),
                 'plot_axis': plot_axis, 'save': True}

    def add_print_variables(self, print_variables, **kwargs):
        for i, name in enumerate(print_variables['names']):
            attr = eval('self.' + name)
            description = string.capwords(name.replace('_', ' '))
            units = print_variables['units'][i]
            sub_names = print_variables.get('sub_names', None)
            if sub_names is not None:
                sub_names = eval(sub_names[i])
            self.add_print_data(attr, description, units=units,
                                sub_names=sub_names, **kwargs)

    @staticmethod
    def combine_print_variables(dict_a, dict_b):
        if dict_b is not None:
            for key in dict_a.keys():
                dict_a[key] += dict_b[key]
        return dict_a

    @classmethod
    def make_name_list(cls):
        name_list = []
        for obj in cls.getinstances():
            obj.name_list = obj.name.split(': ')
            name_list.append(obj.name_list)
        return name_list

    # @classmethod
    # def cluster_objects(cls):
    #     cluster = []


