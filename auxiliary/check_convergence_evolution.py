import os
import json
import timeit
import matplotlib.pyplot as plt
import pemfc
from pemfc.src import simulation
from collections.abc import Iterable as iterable


def build_settings(voltage):
    with open(os.path.join(os.path.dirname(os.path.abspath(pemfc.__file__)), 'settings', 'settings.json')) as file:
        s = json.load(file)

    # Two-Phase and Channel/Land Settings
    s['simulation']['elements'] = 10
    s['simulation']['convergence_criteria'] = 1e-6
    s['simulation']['underrelaxation_factor'] = 0.2
    s['simulation']['maximum_iteration'] = 300
    s['membrane']['underrelaxation_factor'] = 0.2

    s['simulation']['channel_land_discretization'] = False
    s['stack']['calc_temperature'] = True
    s['cathode']['calc_gdl_diffusion'] = False
    s['anode']['calc_gdl_diffusion'] = False
    s['cathode']['calc_two_phase_flow'] = False
    s['anode']['calc_two_phase_flow'] = False

    # Match two-phase source URF to the global thermal URF to prevent runaway
    s['cathode']['two_phase_source_urf'] = s['simulation']['underrelaxation_factor']
    s['anode']['two_phase_source_urf'] = s['simulation']['underrelaxation_factor']

    # Relax the two-phase de-escalation threshold so the scheduler
    # can drop to N=1 updates instead of perpetually shocking the system
    s['cathode']['two_phase_deescalation_error_threshold'] = 5e-4
    s['anode']['two_phase_deescalation_error_threshold'] = 5e-4
    
    # Tune the physical evaporation model to act as an effective macroscopic 
    # mass transfer coefficient for the porous medium, preventing the extreme 
    # latent heat spikes that destabilize the explicit thermal solver.
    for side in ('cathode', 'anode'):
        if 'two_phase_flow' not in s[side]:
            s[side]['two_phase_flow'] = {}
        if 'evaporation_model' not in s[side]['two_phase_flow']:
            s[side]['two_phase_flow']['evaporation_model'] = {}
            
        s[side]['two_phase_flow']['evaporation_model']['type'] = 'HertzKnudsenSchrage'
        s[side]['two_phase_flow']['evaporation_model']['evaporation_coefficient'] = 1e-3
        s[side]['two_phase_flow']['evaporation_model']['condensation_coefficient'] = 1e-3

    s['membrane']['type'] = 'Springer'
    
    s['simulation']['operation_control'] = 'voltage'
    s['simulation']['average_cell_voltage'] = voltage
    
    # Electrochemistry parameters matching the Base Case
    s['cathode']['electrochemistry']['tafel_slope'] = 0.035
    s['cathode']['electrochemistry']['prot_con_cl'] = 3.5
    s['cathode']['electrochemistry']['vol_ex_cd'] = 800000.0
    s['cathode']['electrochemistry']['diff_coeff_cl'] = 1e-7
    s['cathode']['electrochemistry']['diff_coeff_gdl'] = 4e-6
    
    return s

def run_and_plot():
    voltages = [0.5]
    
    fig, axes = plt.subplots(len(voltages), 1, figsize=(8, 10), sharex=True)
    fig.suptitle('Error Evolution at High Current Densities (Two-Phase Flow)', fontsize=14)

    for i, v in enumerate(voltages):
        print(f"\n--- Running Simulation for V = {v} V ---")
        sim = simulation.Simulation(settings=build_settings(v))
        t0 = timeit.default_timer()
        g_data, l_data = sim.run()
        elapsed = timeit.default_timer() - t0
        
        # Extract error history
        curr_errs = l_data[0]['Current Density Error']['value']
        temp_errs = l_data[0]['Temperature Error']['value']
        iterations = list(range(len(curr_errs)))
        
        print(f"Finished V = {v} V | Iters = {len(iterations)} | Time = {elapsed:.1f}s")
        print(f"Final Errors -> Current: {curr_errs[-1]:.2e}, Temp: {temp_errs[-1]:.2e}")

        # Plot on corresponding subplot
        if isinstance(axes, iterable):
            ax = axes[i]
        else:
            ax = axes
        ax.plot(iterations, curr_errs, label='Current Density Error', color='blue', alpha=0.8)
        ax.plot(iterations, temp_errs, label='Temperature Error', color='red', alpha=0.8)
        ax.set_yscale('log')
        ax.set_ylabel('Error')
        ax.set_title(f'Cell Voltage: {v} V')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        ax.legend()
        ax.set_xlabel('Iteration')
    plt.tight_layout()
    
    # Save the plot to the root directory
    plot_path = os.path.join(os.path.dirname(os.path.abspath(pemfc.__file__)), '..', 'auxiliary', 'error_evolution_plot.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to {os.path.abspath(plot_path)}")
    plt.show()

if __name__ == "__main__":
    run_and_plot()