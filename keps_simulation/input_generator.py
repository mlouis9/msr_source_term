# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'common'))

from input_templating import create_inputs
import pandas as pd
import yaml
from string import Template

with open('../source_term/data/simulation_parameters/parameters.yaml', 'r') as f:
    params = yaml.safe_load(f)

rated_power = params['rated_power']
power_history = pd.read_csv('../source_term/data/simulation_parameters/power_history.csv').to_numpy()
powers = power_history[:,1]*rated_power
unique_powers = list(dict.fromkeys(powers))

with open('templates/keps.txt', 'r') as f:
    template = Template(f.read())

for index, power in enumerate(unique_powers):
    templated_str = template.safe_substitute(power=power)
    with open(f'inputs/keps{index}.i', 'w') as f:
        f.write(templated_str)
