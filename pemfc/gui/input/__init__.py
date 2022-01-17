from . import cooling_settings
from . import simulation
from . import cell_settings
from . import manifold_settings
from . import physical_properties
from . import operating_conditions

main_frame_dicts = [cell_settings.tab_dict,
                    manifold_settings.tab_dict,
                    cooling_settings.tab_dict,
                    physical_properties.tab_dict,
                    operating_conditions.tab_dict,
                    simulation.tab_dict]
