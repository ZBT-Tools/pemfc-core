import numpy as np
import pemfc.src.fluid.fluid as fluid
import matplotlib.pyplot as plt

n = 100

fluid_dict = \
    {
        "name": "Cathode Gas Mixture",
        "components": {
            "O2": {
                "state": "gas",
                "molar_fraction": 0.21
            },
            "N2": {
                "state": "gas",
                "molar_fraction": 0.79
            },
            "H2O": {
                "state": "gas-liquid",
                "molar_fraction": 0.0
            }
        },
        "humidity": 0.5,
        "temperature": 303.15,
        "pressure": 101325.0,
        "nodes": n
    }

humid_air = fluid.factory(fluid_dict, backend='pemfc')
temp = np.linspace(0, 400, n)
temp += 273.15
humid_air.update(temperature=temp, pressure=101325)
plt.plot(temp-273.15, humid_air.saturation_pressure/1000.0)
# plt.plot(temp-273.15, humid_air.humidity)
plt.yscale('log')
plt.show()
