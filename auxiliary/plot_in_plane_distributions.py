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
    s['simulation']['elements'] = 25  # Increased for better flow-direction resolution
    s['simulation']['convergence_criteria'] = 1e-6
    s['simulation']['underrelaxation_factor'] = 0.2
    s['simulation']['maximum_iteration'] = 300  
    s['membrane']['underrelaxation_factor'] = 0.2

    s['simulation']['channel_land_discretization'] = True
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
        s[side]['two_phase_flow']['evaporation_model']['evaporation_coefficient'] = 1e-3
        s[side]['two_phase_flow']['evaporation_model']['condensation_coefficient'] = 1e-3

    s['membrane']['type'] = 'Springer'
    s['simulation']['operation_control'] = 'voltage'
    s['simulation']['average_cell_voltage'] = voltage
    
    return s

def run_and_plot():
    voltage = 0.3
    print(f"\n--- Running Simulation for V = {voltage} V ---")
    sim = simulation.Simulation(settings=build_settings(voltage))
    
    t0 = timeit.default_timer()
    sim.run()
    elapsed = timeit.default_timer() - t0
    print(f"Finished V = {voltage} V | Time = {elapsed:.1f}s")
    
    cell = sim.stack.cells[0]
    
    # Extract Saturation at Cathode GDE/MEM interface
    # The 3D array is (nx, ny, nz) -> x=-1 is the GDE/MEM interface
    sat_2d = cell.cathode.two_phase_flow.saturation[-1, :, :]
    
    # Extract Temperature at Cathode GDE/MEM interface from the highly resolved GDL sub-model
    temp_2d = cell.cathode.two_phase_flow.thermal_transport.solution_array[-1, :, :]

    # Extract Current Density (at the MEA) from the global solver
    current_2d = cell.current_density[cell.layer_id['membrane'], :, :]
    
    # Extract Membrane Water Flux from the global solver
    water_flux_2d = cell.membrane.water_flux

    # Physical dimensions for coordinate mapping
    y_len = cell.length # Length in flow direction
    z_len = cell.cathode.flow_field.width_straight_channels / cell.cathode.n_channels # Pitch
    
    # Meshgrids (Note: Thermal grid has nz=2, Two-Phase grid has nz=8)
    Y_sat, Z_sat = np.meshgrid(np.linspace(0, y_len, sat_2d.shape[0]), 
                               np.linspace(0, z_len, sat_2d.shape[1]), indexing='ij')
    Y_tmp, Z_tmp = np.meshgrid(np.linspace(0, y_len, temp_2d.shape[0]), 
                               np.linspace(0, z_len, temp_2d.shape[1]), indexing='ij')
    Y_crs, Z_crs = np.meshgrid(np.linspace(0, y_len, current_2d.shape[0]), 
                               np.linspace(0, z_len, current_2d.shape[1]), indexing='ij')

    fig, axes = plt.subplots(4, 1, figsize=(10, 14), sharex=True)

    c1 = axes[0].pcolormesh(Y_sat, Z_sat, sat_2d, cmap='Blues', shading='auto', vmin=0.0)
    axes[0].set_title('Cathode GDL Liquid Saturation (close to MEM-GDE interface)')
    axes[0].set_ylabel('Width / Pitch [m]')
    fig.colorbar(c1, ax=axes[0], label='Saturation [-]')

    c2 = axes[1].pcolormesh(Y_tmp, Z_tmp, temp_2d, cmap='inferno', shading='auto')
    axes[1].set_title('Temperature at Cathode GDE-MEM Interface')
    axes[1].set_ylabel('Width / Pitch [m]')
    fig.colorbar(c2, ax=axes[1], label='Temperature [K]')

    c3 = axes[2].pcolormesh(Y_crs, Z_crs, current_2d, cmap='viridis', shading='auto')
    axes[2].set_title('Local Current Density')
    axes[2].set_ylabel('Width / Pitch [m]')
    fig.colorbar(c3, ax=axes[2], label='Current Density [A/m²]')

    c4 = axes[3].pcolormesh(Y_crs, Z_crs, water_flux_2d, cmap='coolwarm', shading='auto')
    axes[3].set_title('Membrane Water Flux (Positive = Anode to Cathode)')
    axes[3].set_xlabel('Flow Direction Length [m]')
    axes[3].set_ylabel('Width / Pitch [m]')
    fig.colorbar(c4, ax=axes[3], label='Flux [mol/(s·m²)]')

    plt.tight_layout()
    
    plot_path = os.path.join(os.path.dirname(os.path.abspath(pemfc.__file__)), '..', 'auxiliary', 'in_plane_distributions.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to {os.path.abspath(plot_path)}")
    plt.show()

if __name__ == '__main__':
    run_and_plot()