from string import Template
import pandas as pd
import copy
from pathlib import Path
import periodictable

def element_name(isotope, format_lowercase=False):
    """Get the element name of isotope string in Griffin format"""
    if isotope[-1] == 'M': # Handle metastable isotopes
        isotope = isotope[:-1]
    elif isotope[-2:] == 'M2': # Handle doubly metastable isotopes
        isotope = isotope[:-2]
    translated = isotope.translate(str.maketrans('', '', '0123456789'))
    if format_lowercase and len(translated) == 2:
        translated = translated[0] + translated[1].lower()
    return translated

def lowercase_isotope_name(isotope):
    lowercase_element_name = element_name(isotope, format_lowercase=True)
    return isotope.replace(lowercase_element_name.upper(), lowercase_element_name)

def get_atomic_number(element_symbol):
    """Retrieve the atomic number from the periodictable library."""
    if len(element_symbol) == 2 and element_symbol[-1].isupper():
        element_symbol = element_symbol[0] + element_symbol[1].lower()
    element = getattr(periodictable, element_symbol)
    return element.number  # Atomic number (Z)

def get_zaid(isotope_name):
    # Check for metastable indicators
    metastable = '0'
    if isotope_name[-2:] == 'M2':
        metastable = '2'  # Doubly metastable
        isotope_name = isotope_name[:-2]
    elif isotope_name[-1] == 'M':
        metastable = '1'  # Metastable
        isotope_name = isotope_name[:-1]
    
    # Extract the element symbol and mass number
    element_symbol = element_name(isotope_name)
    mass_number = int(''.join(filter(str.isdigit, isotope_name)))
    
    # Get the atomic number using the periodictable library
    atomic_number = get_atomic_number(element_symbol)
    
    # Format the mass number to always have three digits
    formatted_mass_number = f"{mass_number:03d}"
    
    # Create the ZAID
    zaid = f"{atomic_number}{formatted_mass_number}{metastable}"
    
    return atomic_number, mass_number, zaid

def remove_first_and_last_lines(string):
    lines = string.splitlines()
    lines = lines[1:-1]
    return '\n'.join(lines)

def build_template_string(iterable, template, substitution_dict, top_level=False):
    result = ''
    for index, member in enumerate(reversed(iterable)):
        substitution_dict_copy = copy.deepcopy(substitution_dict)
        for key, expression in substitution_dict_copy.items():
            if key == 'next':
                substitution_dict_copy['next'] = result
            else:
                try:
                    substitution_dict_copy[key] = eval(expression)
                except Exception as e:
                    print(f"The following error occurred when assigning {key} to eval({expression})\n{e}")
        
        result = template.safe_substitute(substitution_dict_copy)

        if not index == len(iterable) - 1:
            result = '\n' + remove_first_and_last_lines(result)
    
    return result

def create_inputs(iterable, template_file_dict, template_dir=Path('templates'), input_dir=Path('inputs'), name_only=True):
    """Create a set of inputs that are templated (for each iterable) according to a template file dict. This dict has keys of template file names, and
    associated substitution dicts that specify what expressions (relating to the iterable) should be substituted for each parameter in the template string."""
    template_files = list(template_file_dict.keys())
    template_files = [ Path(template_file) for template_file in template_files ]
    substitution_dicts = list(template_file_dict.values())
    
    # Read templates into Template objects
    templates = []
    for template_file in template_files:
        with open(template_dir / template_file, 'r') as f:
            template = f.read()
        templates.append(Template(template))


    # Now write template into input files for each iterable
    for template_idx, template_file in enumerate(template_files):
        input_filename = template_file.name.replace(".txt", ".i") if name_only else template_file.replace(".txt", ".i")
        with open(input_dir / input_filename, 'w') as f:
            template = templates[template_idx]
            substitution_dict = substitution_dicts[template_idx]

            text = build_template_string(iterable, template, substitution_dict, top_level=True)
            f.write(text)

def invert_dictionary(original_dict):
    inverted_dict = {}
    
    for key, value_list in original_dict.items():
        for value in value_list:
            inverted_dict[value] = key
    
    return inverted_dict

def process_surrogates(lowercase_isotopes, surrogates_dict):
    elements = set([ element_name(isotope) for isotope in lowercase_isotopes ])
    # Update surrogate dict to include all elements without an explicit surrogate as the Misc element
    surrogate_elements = []
    for element, surrogates in surrogates_dict.items():
        if surrogates == [ 'Misc' ]:
            misc_thermo_element = element
            continue
        surrogate_elements += [ element ] + surrogates
    surrogates_dict[misc_thermo_element] = list( elements - set(surrogate_elements) - {misc_thermo_element} )

    # Get the corresponding surrogate element for each isotope
    inv_surrogate_dict = invert_dictionary(surrogates_dict)
    inv_surrogate_dict.update({element: element for element in surrogates_dict.keys()})
    surrogate_elements = [ element_name(isotope) if element_name(isotope) in surrogates_dict.keys() else inv_surrogate_dict[element_name(isotope)] for isotope in lowercase_isotopes ]

    return surrogates_dict, inv_surrogate_dict, surrogate_elements