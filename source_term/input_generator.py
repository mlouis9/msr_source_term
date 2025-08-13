# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'common'))

from input_templating import create_inputs, get_zaid, element_name, lowercase_isotope_name, process_surrogates
import pandas as pd
from pathlib import Path
from collections import defaultdict
import yaml
import argparse
import json

def read_bateman_output(bateman_csv_file, isotopes, region='core_concentration'):
    data = pd.read_csv(bateman_csv_file)
    zaid_to_isotope = {get_zaid(isotope)[2]: isotope for isotope in isotopes}

    isotopics = {}
    for index, row in data.iterrows():
        if index == 0: # Skip row, contains non-isotope with zaid 3000004
            continue
        try:
            zaid_key = str(int(row['isotope_zaids']))
            isotope = zaid_to_isotope[ zaid_key ]
            isotopics.update({ isotope: row[region] })
        except KeyError: # Non isotope, not tracked by Griffin depletion chain
            pass

    return isotopics

##################################################################
# Constants
##################################################################
PUMP_REMOVAL = 5
input_dir = Path('inputs')
thermofile = 'data/thermo/MSTDB-TC_V4.0_Chlorides_No_Func.dat'
thermo_temp = 923.15 # K
thermo_press = 1.98 # atm

# Initial conditions for Thermochimica (need these to avoid Thermochimica error 5 on first Thermochimica solve)
ics = defaultdict(float)
ics.update({
    'U': 0.33,
    'Cl': 1.66,
    'Na': 0.67,
    'Ar': 2.021,
})

DEFAULT_PARAMS = {
    'dt': 30240,
    'start_time': 0,
    'end_time': 3931200,
    'th_num_steps': 1,
    'initial_isotopics_file': '',
    'power': 150000,
    'index': 0,
}

param_descriptions = {
    'dt': 'Time step size (e.g., 30240)',
    'start_time': 'Simulation start time (e.g., 0)',
    'end_time': 'Simulation end time (e.g., 3931200)',
    'th_num_steps': 'Number of thermal-hydraulic steps (e.g., 5)',
    'initial_isotopics_file': 'The name of the .csv file to read initial isotopics (nuclide densities in 1/(b * cm)). NOTE: This is only intended for reading from the EOL concentrations of a previous depletion simulation',
    'power': 'The thermal power (in W) to run depletion at',
    'index': 'Integer index used to differentiate outputs',
}

##################################################################
# Execution control
##################################################################
parser = argparse.ArgumentParser(
    prog='input_generator',
    description='Generates inputs for source term simulation',
)

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process some control parameters.', formatter_class=argparse.RawTextHelpFormatter)

param_help = "JSON string of control parameters. Expecting the keys:\n"
for key, desc in param_descriptions.items():
    param_help += f"    {key}: {desc}\n"


parser.add_argument('-e', '--external_driver', type=str, help=param_help)

# Parse the arguments
args = parser.parse_args()

# Check if the external_driver argument is provided
if args.external_driver:
    try:
        # Parse the JSON string into a dictionary
        print(args.external_driver)
        control_params = json.loads(args.external_driver)
        print("Control Parameters:", control_params)
        for key, value in control_params.items():
            globals()[key] = value
    except json.JSONDecodeError as e:
        print("Error parsing JSON:", e)
else:
    print("No external driver. Running standalone.")

    # Automatically create variables from dictionary keys
    for key, value in DEFAULT_PARAMS.items():
        globals()[key] = value

if initial_isotopics_file:
    # Read initial isotopics from EOL isotopics of previous simulation
    isotopes = list(pd.read_csv('./data/simulation_parameters/isotopes.csv').keys())
    initial_isotopics = defaultdict(float, read_bateman_output(initial_isotopics_file, isotopes, region='core_concentration'))
    initial_disposal_isotopics = read_bateman_output(initial_isotopics_file, isotopes, region='disposal_concentration')
    initial_disposal_isotopics = '\'' + ' '.join([ f'{key} {value:1.5E}' for key, value in initial_disposal_isotopics.items() ]) + '\''
else:
    initial_isotopics = defaultdict(float, {
        'LI6': 1.28417E-07,
        'LI7': 1.56351E-06,
        'BE9': 1.69193E-06,
        'B10': 3.35001E-08,
        'B11': 1.35692E-07,
        'C12': 2.00777E-05,
        'C13': 2.25364E-07,
        'O16': 1.68790E-05,
        'O17': 6.42932E-09,
        'O18': 3.38385E-08,
        'NA23': 8.43422E-03,
        'MG24': 1.09901E-04,
        'MG25': 1.39133E-05,
        'MG26': 1.53186E-05,
        'AL27': 5.07578E-06,
        'SI28': 3.12093E-07,
        'SI29': 1.58466E-08,
        'SI30': 1.04459E-08,
        'S32': 6.43067E-06,
        'S33': 5.07578E-08,
        'S34': 2.84920E-07,
        'S36': 1.35354E-09,
        'CL35': 1.59544E-02,
        'CL37': 5.10197E-03,
        'K39': 9.46714E-06,
        'K40': 1.18773E-09,
        'K41': 6.83220E-07,
        'CA40': 3.28021E-07,
        'CA42': 2.18929E-09,
        'CA43': 4.56806E-10,
        'CA44': 7.07205E-09,
        'CA46': 1.35350E-11,
        'CA48': 6.32760E-10,
        'TI46': 2.79168E-07,
        'TI47': 2.51759E-07,
        'TI48': 2.49458E-06,
        'TI49': 1.83066E-07,
        'TI50': 1.75283E-07,
        'V50': 4.22981E-09,
        'V51': 1.68770E-06,
        'CR50': 2.94057E-08,
        'CR52': 5.67059E-07,
        'CR53': 6.42999E-08,
        'CR54': 1.60056E-08,
        'MN55': 1.69193E-06,
        'FE54': 9.88931E-08,
        'FE56': 1.55241E-06,
        'FE57': 3.58519E-08,
        'FE58': 4.77123E-09,
        'NI58': 5.75908E-07,
        'NI60': 2.21837E-07,
        'NI61': 9.64398E-09,
        'NI62': 3.07423E-08,
        'NI64': 7.83362E-09,
        'CU63': 2.34061E-07,
        'CU65': 1.04324E-07,
        'ZN64': 1.64557E-06,
        'ZN66': 9.44095E-07,
        'ZN67': 1.38738E-07,
        'ZN68': 6.34472E-07,
        'ZN70': 2.09799E-08,
        'ZR90': 4.35248E-06,
        'ZR91': 9.49170E-07,
        'ZR92': 1.45083E-06,
        'ZR94': 1.47028E-06,
        'ZR96': 2.36870E-07,
        'NB93': 1.69193E-07,
        'MO92': 7.53245E-07,
        'MO94': 4.69509E-07,
        'MO95': 8.08064E-07,
        'MO96': 8.46640E-07,
        'MO97': 4.84737E-07,
        'MO98': 1.22478E-06,
        'MO100': 4.88797E-07,
        'CD106': 8.45963E-10,
        'CD108': 6.02326E-10,
        'CD110': 8.45286E-09,
        'CD111': 8.66266E-09,
        'CD112': 1.63305E-08,
        'CD113': 8.27013E-09,
        'CD114': 1.94436E-08,
        'CD116': 5.06901E-09,
        'SM144': 1.03884E-08,
        'SM147': 5.07239E-08,
        'SM148': 3.80345E-08,
        'SM149': 4.67648E-08,
        'SM150': 2.49728E-08,
        'SM152': 9.05180E-08,
        'SM154': 7.69826E-08,
        'PB204': 9.47478E-09,
        'PB206': 1.63102E-07,
        'PB207': 1.49566E-07,
        'PB208': 3.54628E-07,
        'U234': 4.06631E-05,
        'U235': 3.83891E-03,
        'U236': 1.11971E-05,
        'U238': 2.26123E-04,
    })
    initial_disposal_isotopics = "\'\'"

##################################################################
# Template dicts
##################################################################

template_file_dict_element = {
    'element_aux_kernels.txt': {
        'element_name': 'member[0]',
        'additional_expression': "\' + \' + \'${stability_eps}\' if member[0] == \'Cl\' else \' + \' + \'${ar_moles}\' if member[0] == \'Ar\' else \'\'",
        'next': '',
    },
    'element_aux_variables.txt': {
        'element_name': 'member[0]',
        'next': '',
    },
    'sumset_postprocessors.txt': {
        'element_name': 'member[0]',
        'isotope_sum': 'member[1]',
        'isotope_list': 'member[2]',
        'initial_value': 'member[4]',
        'next': '',
    },

    # Thermochimica template files
    'ics.txt': {
        'element_name': 'member[0]',
        'ic_value': 'member[3][ member[0] ]',
        'next': '',
    },
    'thermo/thermo_transfers.txt': {
        'element_name': 'member[0]',
        'next': '',
    },
    'thermo/thermo_postprocessors.txt': {
        'element_name': 'member[0]',
        'next': '',
    },
}

template_file_dict_isotope = {
    'isotope_postprocessors.txt': {
        'isotope_name': 'member[0]',
        'isotope_zaid': 'member[1]',
        'next': '',
    }, 
    # Total removal functions (Thermochimica + Species transport)
    'total_removal_functions.txt': {
        'isotope_name': 'member[0]',
        'isotope_name_th_removal': 'member[3]',
        'th_removal_name': '\'\' if member[3] == \'1\' else member[3] + \' \'',
        'element_name': 'member[2]',
        'dt': str(dt),
        'next': '',
    },
    'total_removal_postprocessors.txt': {
        'isotope_name': 'member[0]',
        'isotope_name_th_removal': '\'\' if member[3] == \'1\' else member[3] + \' \'',
        'element_name': 'member[2]',
        'next': '',
    }
}

# Template file dict for templates requiring an iterable with a single item
template_file_dict_single = {
    'vector_postprocessors.txt': {
        'initial_isotope_names': 'member[0]',
        'initial_isotope_densities': 'member[1]',
        'isotope_names': 'member[2]',
        'total_removal_rates': 'member[3]',
        'power': str(power),
        'disposal_atomic_densities': initial_disposal_isotopics,
    },
    'thermo/thermo_constants.txt': {
        'thermofile': 'member[4]',
        'thermo_elements_list': 'member[5]',
        'temperature': str(thermo_temp),
        'pressure': str(thermo_press),
    },
    'outputs.txt': {
        'sum_set_list': 'member[6]',
        'gas_densities_list': 'member[9]',
        'thermo_removal_list': 'member[7]',
        'th_removal_list': 'member[8]',
        'index': str(index)
    },
    'th/th_executioner.txt': {
        'th_num_steps': str(th_num_steps),
        'dt': str(dt),
        'system_names': 'member[10]',
    },
    'th/th_functions_single.txt': {
        'th_num_steps': str(th_num_steps),
        'dt': str(dt),
    },
    'th/th_user_objects.txt': {
        'keps_solution': f'\'keps{index}_out.e\'',
    },
    'th/th_mesh.txt': {
        'mesh_file': f'\'keps{index}_out.e\'',
    },
    'dep_executioner.txt': {
        'start_time': str(start_time),
        'end_time': str(end_time),
        'dt': str(dt),
    },
    'dep_controls.txt': {
        'start_time': str(start_time),
        'end_time': str(end_time),
        'index': str(index),
    }
}

template_file_dict_species_isotope = {
    'th/th_transfers.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    },
    'th/th_postprocessors.txt': {
        'isotope_name': 'member[0]',
        'element_name': 'member[1]',
        'next': '',
    },
    'th/th_aux_variables.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    },
    'th/th_aux_kernels.txt': {
        'isotope_name': 'member[0]',
        'element_name': 'member[1]',
        'next': '',
    },
    'th/th_functions.txt': {
        'isotope_name': 'member[0]',
        'element_name': 'member[1]',
        'dt': str(dt),
        'th_num_steps': str(th_num_steps),
        'next': '',
    },
}

template_file_dict_species_element = {
    'th/th_thermo_removal_transfers.txt': {
        'element_name': 'member[0]',
        'next': '',
    },
}

##################################################################
# Write template files
##################################################################

# Read in list of isotopes and elements
thermo_elements = list(pd.read_csv('./data/simulation_parameters/thermo_elements.csv').keys())
with open('./data/simulation_parameters/surrogates.yaml', 'r') as f:
    surrogates_dict = yaml.safe_load(f)

upper_elements = [ element.upper() for element in thermo_elements ]
isotopes = list(pd.read_csv('./data/simulation_parameters/isotopes.csv').keys())
isotopes = [ lowercase_isotope_name(isotope) for isotope in isotopes ]

surrogates_dict, inv_surrogate_dict, surrogate_elements = process_surrogates(isotopes, surrogates_dict)

# -------------------------------
# Create element wise input files
# -------------------------------
element_isotope_list = [ [ isotope for isotope in isotopes if element == element_name(isotope) or element_name(isotope) in surrogates_dict[element] ] \
                        for element in thermo_elements ]
isotope_sums = [ ' + '.join(element_isotopes) for element_isotopes in element_isotope_list ]
isotope_lists = [ ' '.join(element_isotopes) for element_isotopes in element_isotope_list ]
ic_list = [ ics for element in thermo_elements ]
initial_element_densities = [ sum([ initial_isotopics[isotope.upper()] for isotope in element_isotopes ]) for element_isotopes in element_isotope_list ]
iterable = list(zip(thermo_elements, isotope_sums, isotope_lists, ic_list, initial_element_densities))
create_inputs(iterable, template_file_dict_element, input_dir=input_dir)

# -------------------------------
# Create isotope wise input files
# -------------------------------
zaids = [ get_zaid(isotope)[2] for isotope in isotopes ]
species_elements = list(pd.read_csv('./data/simulation_parameters/species_elements.csv').keys())
species_isotopes = [ isotope for isotope in isotopes if element_name(isotope) in species_elements ]
isotope_th_removals = [ f'{isotope}_th_removal' if element_name(isotope) in species_elements else '1' for isotope in isotopes ]
iterable = list(zip(isotopes, zaids, surrogate_elements, isotope_th_removals))
create_inputs(iterable, template_file_dict_isotope, input_dir=input_dir)

# ----------------------------------------------
# Template vectorpostprocessors and executioners
# ----------------------------------------------
initial_isotope_names = ' '.join(list(initial_isotopics.keys())).upper() # Isotope names in bateman postprocessor must be uppercase
initial_isotope_densities = list(initial_isotopics.values())
initial_isotope_densities = ' '.join([ f'{density:1.5E}' for density in initial_isotope_densities ])
isotope_names = ' '.join(isotopes).upper() # Isotope names in bateman postprocessor must be uppercase
total_removal_rates = ' '.join([ f'{isotope}_total_removal_func' for isotope in isotopes ])
gas_densities_list = ' '.join([ f'{isotope}_gas' for isotope in isotopes if element_name(isotope) in species_elements ])
thermo_elements_list = ' '.join(thermo_elements)
sum_set_list = ' '.join([ f'{element}SumSet' for element in thermo_elements ])
thermo_removal_list = ' '.join([ f'{element}_removal' for element in thermo_elements ])
th_removal_list = ' '.join([ f'{isotope}_th_removal' for isotope in isotopes if element_name(isotope) in species_elements ])
species_isotopes_list = ' '.join(species_isotopes)

iterable = [ 
    ( 
        initial_isotope_names, 
        initial_isotope_densities, 
        isotope_names, 
        total_removal_rates,
        thermofile,
        thermo_elements_list,
        sum_set_list,
        thermo_removal_list,
        th_removal_list,
        gas_densities_list,
        species_isotopes_list,
    ) 
]
create_inputs(iterable, template_file_dict_single, input_dir=input_dir)

# --------------------------------------
# Template species transport input files
# --------------------------------------

species_surrogates = [ inv_surrogate_dict[element_name(isotope)] for isotope in species_isotopes ]
iterable = list(zip(species_isotopes,species_surrogates))
create_inputs(iterable, template_file_dict_species_isotope, input_dir=input_dir)

unique_species_surrogates = set(species_surrogates)
iterable = list(zip(unique_species_surrogates,))
create_inputs(iterable, template_file_dict_species_element, input_dir=input_dir)