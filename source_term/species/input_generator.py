# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..', 'common'))

from input_templating import create_inputs, element_name, lowercase_isotope_name, process_surrogates
import pandas as pd
import math
import openmc
import numpy as np
import scipy
import periodictable
from thermo.chemical import Chemical
import yaml

cm2tom2 = 1E-04
cm3tom3 = 1E-06
PastocP = 1E+03
R = 8.314 # J/(mol * K)
TG_FLUXES = np.array([1.92e15, 1.92e15])

def get_boiling_point(element_name):
    chemical = Chemical(element_name)
    return chemical.Tb # Boiling point in Kelvin

def stokes_einstein_diffusion_coef(stokes_radius, dynamic_viscosity, temperature):
    """Computes the diffusion coefficient (in m^2/s) for spherical particles in a low Re liquid given its
    stokes_radius, dynamic_viscosity, and the fluid temperature"""
    kB = 1.380649E-23
    return kB * temperature / ( 6 * math.pi * dynamic_viscosity * stokes_radius )

def wilke_chang_diffusion_coef_gas(association_factor, solvent_M, solvent_mu, T, element, P=101325):
    """Implement the Wilke Chang formula to compute the diffusion coefficient (in m^2/s) for solutes in a liquid solvent
    
    Parameters
    ----------
    association_factor
    solvent_M
        molar mass in g/mol
    solvent_mu
        dynamic viscosity in Pa*s
    T
        temperature in K
    element
        element symbol
    P
        pressure in Pa"""
    Tb = get_boiling_point(element)
    solute_molar_vol_bp = R * Tb/P / cm3tom3 # cm3/mol
    return 7.4E-08 * ( association_factor * solvent_M )**(0.5) * T / ( solvent_mu * PastocP * (solute_molar_vol_bp)**(0.6) ) * cm2tom2

def cumulative_fission_yield(isotope, tg_fluxes, fy_dictionary):
    """Calculates the group averaged fission yields of the specified isotope from fission
    of the specified fissile isotope through the fission yields dictionary, using the given
    two-group fluxes to weight. This assumes the fission yield dictionary has only three
    tabulated energies for each isotope, the first two of which correspond to the thermal
    group, and the latter two to the fast group. This is the case for U-235. It is
    assumed that the thermal flux is at index 0, and the fast flux at index 1. Implements
    the following equation, where $X$ refers to the isotope:
    $$
        \\overline{\\gamma} = \\frac{\\phi_1\\int_{E_1}\\gamma_{X}(E) dE + \\phi_0 \\int_{E_0}\\gamma_X(E)dE}{\\int_{E_1} \\phi_1 dE + \\int_{E_0}\\phi_0 dE}
    $$
    """

    energies = u235.energies
    groups_indices = [ [0, 1], [1, 2] ]
    groups_integral_yield = []
    groups_flux_integral = []
    if isotope[-1] == 'M' or isotope[-2:] == 'M2':
        # Metastable and doubly metastable isotopes are not produced from fisison
        return 0

    for index, group_indices in enumerate(groups_indices):
        phi_g = tg_fluxes[index]
        try:
            group_yield_vals = np.array([ 
                fy_dictionary.cumulative[index][isotope].nominal_value for index in group_indices 
            ])
        except:
            return 0
        group_integral_yield = \
            scipy.integrate.trapezoid(group_yield_vals * phi_g, energies[group_indices])
        groups_integral_yield.append(group_integral_yield)

        group_flux_integral = \
            scipy.integrate.trapezoid(np.repeat(phi_g, len(group_indices)), energies[group_indices])
        groups_flux_integral.append(group_flux_integral)

    return sum(groups_integral_yield) / sum(groups_flux_integral)

# Constants

template_file_dict = {
    'constants.txt': {
        'isotope_name': 'member[0]',
        'molecular_diffusion_coef': 'member[2][ element_name(member[0], format_lowercase=True) ]',
        'fiss_yield': 'member[3]',
        'decay_const': 'member[4]',
        'next': '',
    },
    'variables.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    },
    'fv_kernels.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    },
    'postprocessors.txt': {
        'isotope_name': 'member[0]',
        'element_name': 'member[5]',
        'next': '',
    },
    'functions.txt': {
        'isotope_name': 'member[0]',
        'element_name': 'member[5]',
        'next': '',
    },
    'problem.txt': {
        'isotope_names': 'member[1]',
    },
    # Aux kernels/variables for MultiAppCopyTransfers
    'aux_kernels.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    },
    'aux_variables.txt': {
        'isotope_name': 'member[0]',
        'next': '',
    }
}

template_file_dict_element = {
    'element_postprocessors.txt': {
        'element_name': 'member[0]',
        'next': '',
    }
}

template_file_dict_single = {
    'outputs.txt': {
        'th_removal_list': 'member[0]',
    }
}

# Read in list of isotopes and elements
elements = list(pd.read_csv('../data/simulation_parameters/species_elements.csv').keys())
isotopes = list(pd.read_csv('../data/simulation_parameters/isotopes.csv').keys())
lowercase_isotopes = [lowercase_isotope_name(isotope) for isotope in isotopes]
# Read in and process surrogate data
thermo_elements = list(pd.read_csv('../data/simulation_parameters/thermo_elements.csv').keys())
with open('../data/simulation_parameters/surrogates.yaml', 'r') as f:
    surrogates_dict = yaml.safe_load(f)

_, inv_surrogate_dict, _ = process_surrogates(lowercase_isotopes, surrogates_dict)

isotopes = [ isotope for isotope in lowercase_isotopes if element_name(isotope) in elements ]
isotope_list = [ ' '.join(isotopes) for isotope in isotopes ]
print(f"Number of isotopes: {len(isotopes)}")

# Stokes radii for each element (in m)
stokes_radii = {
    'Kr': 180E-12, # Obtained from kinetic diameter: https://en.wikipedia.org/wiki/Kinetic_diameter
    'Xe': 198E-12, # Obtained from kinetic diameter: https://en.wikipedia.org/wiki/Kinetic_diameter
    'Cs': 343E-12, # Obtained from Van Der Waals Radius: https://en.wikipedia.org/wiki/Caesium
    'Cl': 343E-12,
}

# ----------------------------------------------------------------------------------------
# Salt properties necessary for Wilke-Chang estimation of molecular diffusion coefficients
# ----------------------------------------------------------------------------------------
association_factor = 2.26 # Conservatively estimate as that for water

# Masses (g/mol), obtained from JANIS
M_Na   = 22.99
M_Cl   = 35.453
M_Mg   = 24.305
M_U232 = 232.037146
M_U234 = 234.040946	
M_U235 = 235.043923	
M_U236 = 236.045562
M_U238 = 238.050783

# Composition (atomic) fractions of fuel salt, obtained from: https://collab.terrapower.com/display/MCREX/RCS-CALC-0014+%7C%7C+Molecular+Composition+of+Fuel+Salt+Diluted+with+Flush+Salt Table 9-1
af_NaCl = 0.6648
af_UCl3 = 0.3245
af_MgCl2 = 0.0107 # flush salt

# Composition (weight) fractions of Uranium isotopes, obtained from: https://collab.terrapower.com/display/MCREX/FUEL-EDD-0001+%7C%7C+NaCl-UCl3+Eutectic+Salt+Thermophysical+Properties Table 8-1
wf_U232 = 1.76E-10
wf_U234 = 9.83E-03
wf_U235 = 9.32E-01
wf_U236 = 2.73E-03
wf_U238 = 5.56E-02

# Convert the weight fractions to atomic fractions via the expression $af_i = \frac{\frac{wf_i}{M_i}}{\sum_i \frac{wf_i}{M_i}}$
denom   = wf_U232/M_U232 + wf_U234/M_U234 + wf_U235/M_U235 + wf_U236/M_U236 + wf_U238/M_U238
af_U232 = ( wf_U232 / M_U232 ) / denom
af_U234 = ( wf_U234 / M_U234 ) / denom
af_U235 = ( wf_U235 / M_U235 ) / denom
af_U236 = ( wf_U236 / M_U236 ) / denom
af_U238 = ( wf_U238 / M_U238 ) / denom

M_U = af_U232 * M_U232 + af_U234 * M_U234 + af_U235 * M_U235 + af_U236 * M_U236 + af_U238 * M_U238
M_tot = af_NaCl * (M_Na + M_Cl) + af_UCl3 * (M_U + 3*M_Cl) + af_MgCl2 * (M_Mg + 2 * M_Cl)
# ----------------------------------------------------------------------------------------

# Molecular diffusion coefficients for each element computed from the Stokes-Einstein relation
dynamic_viscosity = 4.9E-3 # Pa*s Obtained from https://collab.terrapower.com/display/MCREX/FUEL-EDD-0001+%7C%7C+NaCl-UCl3+Eutectic+Salt+Thermophysical+Properties, Table 9-1 (L_Viscosity_650)
T = 923.15
diffusion_coefficients = { element: wilke_chang_diffusion_coef_gas(association_factor, M_tot, dynamic_viscosity, T, element) for element in elements }
# diffusion_coefficients = { element: stokes_einstein_diffusion_coef(stokes_radius, dynamic_viscosity, T) for element, stokes_radius in stokes_radii.items() }
diffusion_coefficients = [ diffusion_coefficients for isotope in isotopes ]

# Calculate (cumulative) fission yields for each isotope
u235 = openmc.data.FissionProductYields('../data/dt/nfy-092_U_235.endf')
fission_yields = [ cumulative_fission_yield(lowercase_isotope_name(isotope), TG_FLUXES, u235) for isotope in isotopes ]

# Get decay constants for each isotope
decay_constants = [ openmc.data.decay_constant(lowercase_isotope_name(isotope)) for isotope in isotopes ]
species_surrogates = [ inv_surrogate_dict[element_name(isotope)] for isotope in isotopes ]

iterable = list(zip(
    isotopes, 
    isotope_list, 
    diffusion_coefficients, 
    fission_yields, 
    decay_constants,
    species_surrogates,
))
create_inputs(iterable, template_file_dict)

#############################
# Create Element-Wise Inputs
#############################

species_surrogates_set = list(set(species_surrogates))
iterable = list(zip(species_surrogates_set))
create_inputs(iterable, template_file_dict_element)

#######################
# Create Single Inputs
#######################

th_removal_list = ' '.join([ isotope + '_th_removal' for isotope in isotopes ])
iterable = [
    (
        th_removal_list,
    )
]
create_inputs(iterable, template_file_dict_single)