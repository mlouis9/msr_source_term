# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../', 'common'))

from input_templating import lowercase_isotope_name

import openmc
import openmc.deplete
import openmc.mgxs as mgxs
import subprocess
import shutil
import numpy as np
import pandas as pd
import pickle
import yaml
import itertools

# -----------------------------------------------------
# Get two group flux spectrum from depletion simulation
# -----------------------------------------------------
# take at initial timestep since spectrum doesn't change much
statepoint = openmc.StatePoint('openmc_simulation_n0.h5')
tally = statepoint.get_tally(id=95)
np.savetxt('tg_fluxes.txt', tally.mean.flatten())
shutil.copy('tg_fluxes.txt', '../source_term/data/simulation_parameters/tg_fluxes.txt')

# -----------------------
# Generate cross sections
# -----------------------
griffin_xs_gen_command = "ml use.moose griffin-openmpi; griffin-opt --isoxml-input mgxs_1g_control.xml"
subprocess.run(griffin_xs_gen_command, shell=True)
shutil.copy('mgxs.xml', '../source_term/data/simulation_parameters/mgxs_1g.xml')

# --------------------------
# Get fission rate densities
# --------------------------
# Load params
with open('../source_term/data/simulation_parameters/parameters.yaml', 'r') as f:
    params = yaml.safe_load(f)

joule_per_ev = 1.60218e-19
# Scale power to the fraction of the total fuel volume that was explicitly modeled
power_history = pd.read_csv('../source_term/data/simulation_parameters/power_history.csv')
power_history = power_history.to_numpy()
rated_power = params['rated_power']
times = power_history[:,0]
powers = rated_power * power_history[:, 1]

statepoint_files = [f'openmc_simulation_n{i}.h5' for i in range(len(times))]

# Retrieve geometry and library objects
geometry = openmc.Geometry.from_xml('geometry.xml')
with open('mgxs/mgxs.pkl', 'rb') as f:
    library = pickle.load(f)

# Load fuel cell object
material_cells = geometry.get_all_material_cells()
fuel_cell_id = [ key for key, cell_object in material_cells.items() if cell_object.name == 'reactor_cell' ][0]
fuel_cell = material_cells[fuel_cell_id]

# Read materials and get available isotopes
materials = openmc.Materials.from_xml('materials.xml')
fuel = [mat for mat in materials if mat.name == 'HEUF'][0]
isotopes = list(pd.read_csv('../source_term/data/simulation_parameters/isotopes.csv').keys())
data_library = openmc.data.DataLibrary.from_xml()
data_library_nuclides = [ entry['materials'][0] for entry in data_library ]
isotopes = library.nuclides

fuel_vol = fuel.volume

cross_sections = { isotope: [] for isotope in isotopes }
fiss_rate_densities = []
for depletion_index, statepoint_file in enumerate(statepoint_files):
    # Load the StatePoint file
    sp = openmc.StatePoint(statepoint_file)
    library.load_from_statepoint(sp)
    
    # Get heating tally
    heating_tally = sp.get_tally(name='heating')
    heating_tally_df = heating_tally.get_pandas_dataframe()
    heating_rate_ev = heating_tally.mean.flatten()[0]

    # Get fission rate tally
    fission_rate_tally = sp.get_tally(name='fission_rate')
    fission_rate_tally_df = fission_rate_tally.get_pandas_dataframe()
    fission_rate = fission_rate_tally.mean.flatten()[0]

    # Compute the source rate
    source_per_sec = powers[depletion_index] / ( heating_rate_ev * joule_per_ev )
    fission_rate_sec = fission_rate * source_per_sec
    fiss_rate_densities.append(fission_rate_sec/fuel_vol)

data = {'Timesteps': np.array(times), 'Fission Rate Density': np.array(fiss_rate_densities)}
df = pd.DataFrame(data=data)
df.to_csv('fiss_rate_densities.csv', index=False)
shutil.copy('fiss_rate_densities.csv', '../source_term/data/simulation_parameters')