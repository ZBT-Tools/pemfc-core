import numpy as np
import timeit
import src.simulation as simulation
import sys
import os

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


def main():
    np.seterr(all='raise')
    start_time = timeit.default_timer()
    sim = simulation.Simulation()
    sim.timing['start'] = start_time
    sim.timing['initialization'] = timeit.default_timer()
    # simulation.timing['start'] = start_time
    g_data, l_data = sim.run()
    # sim.output.print_global_data(sim, g_data)
    sim.output.save_settings()
    return g_data, l_data, sim


if __name__ == "__main__":
    global_data, local_data, sim = main()

