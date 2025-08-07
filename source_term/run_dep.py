import subprocess
import numpy as np
from string import Template
from math import ceil
import argparse
import pandas as pd
import yaml

parser = argparse.ArgumentParser(
    prog='run_bateman.py',
    description='Runs a sequence of Bateman depletion calculations at different power levels (automatically normalized)',
    epilog=''
)
parser.add_argument('-p', '--procs', default=1, type=int)
parser.add_argument('-c', '--cores', default=1, type=int)
parser.add_argument('-io', '--input_only', default=-1, type=int)
args = parser.parse_args()

num_procs = args.procs
num_cores = args.cores
input_only = args.input_only
bluecrab_run_command = f"""
ml use.moose moose-dev-openmpi/2025.07.22
mpiexec -n {num_procs} moose-dev-exec "/projects/MCRE_studies/louime/bin/blue_crab-opt_v4 -i run_dep.i --n-threads={num_cores}"
"""

input_generator_run_command_template = Template(
"""python input_generator.py -e '{"dt": ${dt}, "start_time": ${start_time}, """
""""end_time": ${end_time}, "th_num_steps": ${th_num_steps}, """
""""initial_isotopics_file": \"${initial_isotopics_file}\", "power": ${power}, "index": ${index}}'"""
)

species_run_command = "cd species; python input_generator.py; cd .."

def run_bluecrab():
    subprocess.run(bluecrab_run_command, shell=True)

def generate_species_inputs():
    subprocess.run(species_run_command, shell=True)

def generate_inputs(dt, start_time, end_time, th_num_steps, initial_isotopics_file, power, index):
    input_generator_run_command = input_generator_run_command_template.safe_substitute({
        'dt': dt,
        'start_time': start_time,
        'end_time': end_time,
        'th_num_steps': th_num_steps,
        'initial_isotopics_file': initial_isotopics_file,
        'power': power,
        'index': index,
    })
    subprocess.run(input_generator_run_command, shell=True)

with open('data/simulation_parameters/parameters.yaml', 'r') as f:
    params = yaml.safe_load(f)
rated_power = params['rated_power']
dts = params['dts']
th_num_steps = params['th_num_steps']

power_history = pd.read_csv('data/simulation_parameters/power_history.csv')
power_history = power_history.to_numpy()

# Create timesteps.csv
timesteps = power_history[:,0]
timesteps = np.concatenate([ np.arange(timesteps[i], timesteps[i+1], dts[i]) for i in range(len(timesteps) - 1) ] + [np.array([timesteps[-1]])])
midpoints = np.concatenate([np.array([0]), (timesteps[1:] + timesteps[:-1])/2])
df = pd.DataFrame.from_dict({'Timesteps': timesteps, 'Depletion Values': midpoints})
df.to_csv('timesteps.csv', index=False)

power_history = power_history[1:, :]
generate_species_inputs()
for index, row in enumerate(power_history):
    if index == 0:
        start_time = 0
        initial_isotopics_file = ''
    else:
        start_time = power_history[index-1, 0]
        previous_start_time = power_history[index-2,0] if index >= 2 else 0

        previous_delta_t = start_time - previous_start_time
        num_timesteps = ceil(previous_delta_t/dts[index-1])
        initial_isotopics_file = f'results/isotopics_{index-1}_bateman_{num_timesteps:04d}.csv'
    end_time = power_history[index, 0]
    dt = dts[index]
    power = rated_power * row[1]
    if input_only != -1:
        if index == input_only:
            generate_inputs(dt, start_time, end_time, th_num_steps, initial_isotopics_file, power, index)
    else:
        generate_inputs(dt, start_time, end_time, th_num_steps, initial_isotopics_file, power, index)
        run_bluecrab()