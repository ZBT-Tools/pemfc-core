import os
import json
from pemfc.main_app import main

# load settings
with open(os.path.join('tests', 'settings.json')) as file:
    ref_settings = json.load(file)

# load settings
with open(os.path.join('pemfc', 'settings', 'settings.json')) as file:
    main_settings = json.load(file)

# load reference results
with open(os.path.join('tests', 'results.json')) as file:
    ref_results = json.load(file)


def test_inputs():
    assert main_settings == ref_settings


def test_global_results():
    g_data, l_data, sim = main(ref_settings)
    results = g_data[0]
    test_result = round(results['Stack Voltage']['value'], 3)
    ref_result = round(ref_results['Stack Voltage']['value'], 3)
    assert test_result == ref_result
