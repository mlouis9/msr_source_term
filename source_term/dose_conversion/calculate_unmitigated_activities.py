# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..', 'modules'))

from input_templating import lowercase_isotope_name, element_name, process_surrogates, get_zaid
import pandas as pd
import numpy as np
from math import ceil
from collections import defaultdict
import openmc
import yaml
import matplotlib.pyplot as plt

MAX_REMOVAL_RATE = 1
perbcm_to_percm3 = 1E+24
Bq_to_Ci = 3.7E+10

power_history = pd.read_csv('../power_history.csv').to_numpy()[1:, :]
N_pow_steps = power_history.shape[0]
dts = [32400, 34200, 3600]

# Process surrogates
isotopes = list(pd.read_csv('../data/simulation_parameters/isotopes.csv').keys())
isotopes = [ lowercase_isotope_name(isotope) for isotope in isotopes ]
with open('../data/simulation_parameters/surrogates.yaml', 'r') as f:
    surrogates_dict = yaml.safe_load(f)

surrogates_dict, inv_surrogate_dict, surrogate_elements = process_surrogates(isotopes, surrogates_dict)

# First read removal fractions and calculate removal deay constants
lambda_df_list = []
ar_removal_list = []
for index in range(N_pow_steps):
    # Read removal fractions from .csv
    thermo_removal_file = f'../run_dep_thermo_removal_fractions_{index}.csv'
    th_removal_file = f'../run_dep_th_removal_fractions_{index}.csv'
    thermo_df = pd.read_csv(thermo_removal_file, index_col=False)
    th_df = pd.read_csv(th_removal_file, index_col=False)

    # Process removal fractions into removal decay constants
    zero_data = np.zeros_like(th_df.to_numpy(dtype=np.float64))
    times_col = th_df.to_numpy(dtype=np.float64)[:, 0]
    zero_data[:, 0] = times_col
    lambda_df = pd.DataFrame(zero_data, columns=['time'] + isotopes)
    lambda_df['time'] = times_col
    ar_removal_list.append(thermo_df['Ar_removal_post'])
    for isotope, surrogate in zip(isotopes, surrogate_elements):
        th_removal = th_df[f'{isotope}_th_removal_post']
        thermo_removal = thermo_df[f'{surrogate}_removal_post']
        both_unity = np.logical_and(th_removal.to_numpy() == 1, thermo_removal.to_numpy() == 1)
        lambda_df.loc[both_unity, isotope] = MAX_REMOVAL_RATE
        lambda_df.loc[~both_unity, isotope] = -np.log(1 - th_removal[~both_unity] * thermo_removal[~both_unity]) / dts[index]

    if index == 0:
        lambda_df_list.append(lambda_df)
    else:
        lambda_df_list.append(lambda_df[1:]) # Skip initial time (zero values, initialization)

lambda_df = pd.concat(lambda_df_list, ignore_index=True)

# Read isotpics at each timestep
isotopics_df_list = []
timestep_dts = []
for index in range(N_pow_steps):
    start_time = 0 if index == 0 else power_history[index-1, 0]
    end_time = power_history[index,0]
    delta_t = end_time - start_time
    num_timesteps = ceil(delta_t/dts[index])
    for time_index in range(num_timesteps+1):
        if time_index == 0 and not index == 0: # Skip beginning of timestep isotopics files (redundant from previous simulation)
            continue
        timestep_dts.append(dts[index]) # Store dt's for each timestep
        isotopics_file = f'../results/isotopics_{index}_bateman_{time_index:04d}.csv'
        isotopics_df = pd.read_csv(isotopics_file, index_col=False)
        new_isotopics_df = pd.DataFrame(np.zeros((1, len(isotopes) + 1)), columns=['time'] + isotopes)
        new_isotopics_df['time'] = start_time + dts[index] * time_index
        for isotope in isotopes:
            new_isotopics_df[isotope] = isotopics_df[isotopics_df['isotope_zaids'] == int(get_zaid(isotope)[2])]['core_concentration'].to_numpy().flatten()
        isotopics_df_list.append(new_isotopics_df)

timestep_dts = np.array(timestep_dts)
isotopics_df = pd.concat(isotopics_df_list, ignore_index=True)

# Compute upper bound on (undecayed) removal
off_gas_isotopics = defaultdict(float)
for isotope in isotopes:
    isotopics = isotopics_df[isotope].to_numpy()
    max_N = np.maximum(isotopics[:-1], isotopics[1:]) * perbcm_to_percm3 # Conver to 1/cm3
    lambdas = lambda_df[isotope].to_numpy()[:-1] # Use previous step values for lambda, based on execution precedence of simulation. Actually used in bateman
    off_gas_isotopics[isotope] += np.sum(lambdas * max_N * timestep_dts[1:])

# Convert isotopics to activities (in Ci)
off_gas_activities = {
    isotope: [ openmc.data.decay_constant(lowercase_isotope_name(isotope)) * off_gas_isotopics[isotope] / Bq_to_Ci ] 
    for isotope in isotopes 
}
df = pd.DataFrame.from_dict(off_gas_activities)
df.to_csv('off_gas_activities_unmitigated.csv', index=False)

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
plt.savefig('top_n_isotopes_activity_unmitigated.png', dpi=500)