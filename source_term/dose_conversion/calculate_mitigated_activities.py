# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..', 'modules'))

from input_templating import get_zaid, lowercase_isotope_name
import pandas as pd
import numpy as np
from math import ceil
from collections import defaultdict
import openmc
import matplotlib.pyplot as plt

perbcm_to_percm3 = 1E+24
Bq_to_Ci = 3.7E+10
power_history = pd.read_csv('../power_history.csv').to_numpy()[1:, :]
isotopes = list(pd.read_csv('../data/simulation_parameters/isotopes.csv').keys())
N_pow_steps = power_history.shape[0]
dts = [32400, 34200, 3600]

df_list = []
off_gas_isotopics = defaultdict(float)
index = N_pow_steps - 1 # Read from last isotopics file
start_time = power_history[index-1, 0] if index != 0 else 0
end_time = power_history[index, 0]
num_timesteps = int((end_time - start_time)/dts[index])

isotopics_file = f'../results/isotopics_{index}_bateman_{num_timesteps:04d}.csv'
df = pd.read_csv(isotopics_file, index_col=False)
for isotope in isotopes:
    isotope_data = df[df['isotope_zaids'] == int(get_zaid(isotope)[2])] * perbcm_to_percm3 # Convert to 1/cm3
    disposal_concentration = isotope_data['disposal_concentration'].to_numpy().flatten()[0]
    off_gas_isotopics[isotope] += disposal_concentration

# Convert isotopics to activities
off_gas_activities = {
    isotope: [ openmc.data.decay_constant(lowercase_isotope_name(isotope)) * off_gas_isotopics[isotope] / Bq_to_Ci ] 
    for isotope in isotopes 
}
df = pd.DataFrame.from_dict(off_gas_activities)
df.to_csv('off_gas_activities_mitigated.csv', index=False)

# Also create a histogram plot of top n isotopes
top_n = 20
top_n_isotopes = df.T.nlargest(top_n, 0)
plt.figure(figsize=(12, 8))
plt.bar(top_n_isotopes.index, top_n_isotopes[0])
plt.xlabel('Isotope')
plt.ylabel('Activity (Ci)')
plt.title(f'Top {top_n} Isotopes by Activity')
plt.xticks(rotation=45, ha='right')
plt.grid()
plt.yscale('log')
plt.tight_layout()
plt.savefig('top_n_isotopes_activity_mitigated.png', dpi=500)
