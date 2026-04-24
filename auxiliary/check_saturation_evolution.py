import os
import json
import timeit
import numpy as np
import matplotlib.pyplot as plt
import pemfc
from pemfc.src import simulation

def build_settings(voltage):
    with open(os.path.join(os.path.dirname(os.path.abspath(pemfc.__file__)), 'settings', 'settings.json')) as file:
        s = json.load(file)

    # Two-Phase and Channel/Land Settings
    s['simulation']['channel_land_discretization'] = True
    s['simulation']['elements'] = 10
    s['simulation']['convergence_criteria'] = 1e-6
    s['simulation']['underrelaxation_factor'] = 0.05
    s['simulation']['maximum_iteration'] = 300
    s['stack']['calc_temperature'] = True
    
    s['cathode']['calc_gdl_diffusion'] = True
    s['anode']['calc_gdl_diffusion'] = True
    s['cathode']['calc_two_phase_flow'] = True
    s['anode']['calc_two_phase_flow'] = True

    # Match two-phase source URF to the global thermal URF to prevent runaway
    s['cathode']['two_phase_source_urf'] = s['simulation']['underrelaxation_factor']
    s['anode']['two_phase_source_urf'] = s['simulation']['underrelaxation_factor']

    # Relax the two-phase de-escalation threshold
    s['cathode']['two_phase_deescalation_error_threshold'] = 5e-4
    s['anode']['two_phase_deescalation_error_threshold'] = 5e-4
    
    # Tuned macroscopic evaporation/condensation parameters
    for side in ('cathode', 'anode'):
        if 'two_phase_flow' not in s[side]:
            s[side]['two_phase_flow'] = {}
        if 'evaporation_model' not in s[side]['two_phase_flow']:
            s[side]['two_phase_flow']['evaporation_model'] = {}
            
        s[side]['two_phase_flow']['evaporation_model']['type'] = 'HertzKnudsenSchrage'
        s[side]['two_phase_flow']['evaporation_model']['evaporation_coefficient'] = 1e-2
        s[side]['two_phase_flow']['evaporation_model']['condensation_coefficient'] = 1e-2

    s['membrane']['type'] = 'Springer'
    s['simulation']['operation_control'] = 'voltage'
    s['simulation']['average_cell_voltage'] = voltage
    
    return s

def run_and_plot():
    voltages = [0.3, 0.2, 0.1]
    
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    fig.suptitle('Cathode GDL Saturation Evolution over Iterations', fontsize=14)

    for i, v in enumerate(voltages):
        print(f"\n--- Running Simulation for V = {v} V ---")
        sim = simulation.Simulation(settings=build_settings(v))
        
        sat_history = []
        
        # Monkey-patch the stack update to extract saturation history cleanly
        original_update = sim.stack.update
        def patched_update(*args, **kwargs):
            original_update(*args, **kwargs)
            # Extract average and maximum saturation from the Cathode GDL
            cathode_sat = sim.stack.cells[0].cathode.saturation
            sat_history.append((np.average(cathode_sat), np.max(cathode_sat)))
            
        sim.stack.update = patched_update
        
        t0 = timeit.default_timer()
        g_data, l_data = sim.run()
        elapsed = timeit.default_timer() - t0
        
        avg_sat_hist = [s[0] for s in sat_history]
        max_sat_hist = [s[1] for s in sat_history]
        iterations = list(range(len(sat_history)))
        
        print(f"Finished V = {v} V | Iters = {len(iterations)} | Time = {elapsed:.1f}s")
        print(f"Final Saturation -> Avg: {avg_sat_hist[-1]:.4f}, Max: {max_sat_hist[-1]:.4f}")

        # Plot on corresponding subplot
        ax = axes[i]
        ax.plot(iterations, avg_sat_hist, label='Average Saturation', color='blue', alpha=0.8)
        ax.plot(iterations, max_sat_hist, label='Maximum Saturation', color='orange', linestyle='--', alpha=0.8)
        ax.set_ylabel('Liquid Saturation [-]')
        ax.set_title(f'Cell Voltage: {v} V')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        ax.legend()

    axes[-1].set_xlabel('Iteration')
    plt.tight_layout()
    
    # Save the plot
    plot_path = os.path.join(os.path.dirname(os.path.abspath(pemfc.__file__)),  '..', 'auxiliary', 'cathode_gdl_saturation_evolution_plot.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to {os.path.abspath(plot_path)}")
    plt.show()

if __name__ == "__main__":
    run_and_plot()