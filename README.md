# Introduction
This project demonstrates an enegineering-scale multiphysics technique for Molten Salt Reactor source term calculation using the MOOSE framework. This technique is applied to an open source model of the MCRE experiment at INL, the LOTUS Molten Chloride Reactor (LMCR), detailed [here](https://mooseframework.inl.gov/virtual_test_bed/msr/lotus/lotus_description.html). The calculation sequence integrates thermal hydrualics, species transport, and depletion to calculate the total amount of fission products vented into the off gas system over a typical operating period.

# Directory description
1. [keps_simulation](keps_simulation): Contains the thermal hydraulics model necessary for generating the flow field that is used in species transport calculations. This consists of a $k$-$\epsilon$ MOOSE Navier Stokes simulation of the entire primary loop.
2. [source_term](source_term): Contains the primary source term calculation along with some necessary postprocessing scripts. This consists of a coupled depletion-thermochemistry-species transport calculation that uses the flow field computed from the `keps_simulation` and the multigroup cross sections computed from `xs_gen`. 
3. [xs_gen](xs_gen): Contains the open-source OpenMC neutronics model of the full primarly loop of the LMCR necessary for generating 1 group cross sections for depletion calculations in the main multiphysics depletion simulation. Note there are also tools for converting the cross sections from `.h5` to Griffin's ISOXML format that require either a bluecrab or griffin executable.
4. [common](common): Common python code used in templating input files for the multiphysics simulation. As it stands, the syntax of the experimental depletion feature requires very verbose input, so this is necessary for the time being. You can read more about this in the [common](common) directory.

# Requirements
- **BLUECRAB**: Access to a BLUECRAB executable, which is export controlled and requires a formal request. More information can be found [here]([here](https://inl.gov/ncrc/how-ncrc-works/)).
- **Griffin Depletion with Isotope Removal**: In this project an experimental version of Griffin with a depletion solver that allows for isotope removal is used. Specifically, a `BatemanMultiVPP` vector postprocessor is used, and requires a BLUECRAB executable to be built with this experimental Griffin feature. This feature will eventually be added (though likely not with the same input syntax) to a future version of Griffin.
- **Thermochimica**: Finally, note that your BlueCRAB executable must be built with Thermochimica. Instructions on how to do this are given [here](https://mooseframework.inl.gov/source/actions/ChemicalCompositionAction.html).
- **Python libraries**:
  - pandas
  - periodictable
  - pyyaml
  - thermo
- **OpenMC**: A working installation of OpenMC (with the nuclear data libraries) is required for generating the cross sections in `xs_gen`, guidance is given [here](https://docs.openmc.org/en/stable/quickinstall.html).
- The latest version of **MSTDB-TC** (Molten Salt Thermal Properties Database - Thermo Chemical). Can be obtained [here]([here](https://mstdb.ornl.gov/)). Then place the FactSage `.dat` file for the Chloride database (e.g. `MSTDB-TC_V4.0_Chlorides_No_Func.dat`) in `source_term/data/thermo`.

# Input Templating
This primarily serves to describe the logic of the input templating implemented by `common/intput_templating.py`. Several of the simulations are either required to be run multiple times, or with very verbose, possibly flexible input, and to cope with this a simple python templating system was devised. For the `source_term` simulation, each app has a main input file (and so does each subapp) and they include files that are generated programmatically from the information in the `templates` directory and placed in `inputs`. In `keps_simulation`, the main input files are generated in their entirety by the input templater. The `templates` `.txt` files are written in the form of python Template strings, where `${variable}` indicates a named python variable that will be templated via substituting some dictionary `{'variable': 0.5}`. When a particular pattern is repeated in the input multiple times, e.g.
```
[AuxKernels]
  [Xe135_gas]
    type = ParsedAux
    expression = 'Xe135 * Xe_removal'
    functor_names = 'Xe135 Xe_removal'
    functor_symbols = 'Xe135 Xe_removal'
    variable = Xe135_gas
  []
  [Xe136_gas]
    type = ParsedAux
    expression = 'Xe136 * Xe_removal'
    functor_names = 'Xe136 Xe_removal'
    functor_symbols = 'Xe136 Xe_removal'
    variable = Xe136_gas
  []
```
etc. it is shorted by using
```
[AuxKernels]
  [${isotope_name}_gas]
    type = ParsedAux
    expression = '${isotope_name} * ${element_name}_removal'
    functor_names = '${isotope_name} ${element_name}_removal'
    functor_symbols = '${isotope_name} ${element_name}_removal'
    variable = ${isotope_name}_gas
  []${next}
[]
```
where the keyword `${next}` denotes the next 'block' in the input section, which will be templated with different values of the variables `isotope_name`, `element_name`, as in `source_term/templates/th/th_aux_kernels.txt`. In the main simulation, input templating is controlled by `input_generator.py` and reads data from the `data/simulation_parameters` directory primarily to template the input files for a given set of chemical surrogates (for thermochemistry calculations), power history, etc. If you need to change any of the inputs for your work, you must make edits in three locations
1. The template `.txt` file
2. The main simulation file (change the include if making a new template or changing its name)
3. The `input_generator.py` file if you add any new variables, this script needs to know what to template for those values.
In the input files are generally templated by the type of data that they contain, so, sections (e.g. AuxKernels) that are repeated for every isotope are grouped, those that are repeated for every element are grouped separately, and those that are repeated only once are also separate. For each class, there is a template dictionary created that describes how each of the files in that class should be tempalted in terms of some common data that is given to the `input_templating create_inputs` function. For example, the following file is templated for every element
```
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
```
The keyword `member[i]` refers to the ith member of a tuple that is common among all template files in this class. For the example above, each of the files would be templated for each element via a common list of tuples, where the tuples correspond to information, e.g. the element name, for each element.

# Simulation Parameters
The physics of the source term simulation is parameterized by a few things (each contained in `source_term/data/simulation_parameters`)
1. `power_history.csv`: Specifies the power history of the reactor, in a particular format that is compatible with how MOOSE reads `.csv` files. Note that
    ```
    Timesteps,Power Factor
    0,0.006666666666666667
    3240000,0.006666666666666667
    3924000,0.13333333333333333
    3931200,1.0
    ```
    refers to a power fraction of 0.006666666666666667 of rated power until 3240000 seconds, and a power of 13333333333333333 until 3924000 seconds, and full power from there onward. The rated power is set in `sorce_term/data/simulation_parameters/parameters.yaml` with the `rated_power` parameter. This is currently set to the rated power of the LMCR model.
2. `isotopes.csv`: The list of isotopes to be tracked with explicit removal rates in the depletion solver. By default, all isotopes in the Griffin depletion chain will be tracked by default (and cross sections for each isotope should be generated properly by the `xs_gen` step), but only those specified here will have nonzero removal rates from the salt. In addition, these isotopes will not be tracked for the purposes of thermochemsitry (will not be assigned any element or surrogate).
3. `species_elements.csv`: The list of elements to be included in species transport calculations. This will currently include all of the isotopes of this elemenet specified in `isotopes.csv`. In the future it may be more efficient to group isotopes by half-lives and molecular diffusivity. For conservativism (in estimating the source term), any element not explicitly modeled will have a thermal hydraulic removal factor of 1.
4. `thermo_elements.csv`: The elements to be included in thermochemistry calculations. Any isotopes of elements not included in this list will be assumed to be an isotope of an element in this list by the mapping defined in the `surrogates.yaml` file. Note that this must include at least the Argon cover gas in the pump, and the elements must all be contained in the MSTDB-TC `.dat` file you're using, or Thermochimica will give an error.
5. `surrogates.yaml`: This is a mapping (a dictionary in `.yaml` syntax) keyed by thermochemical elements with values equal to the elements that will be mapped to this element (key) as a surrogate. For example
    ```
    U:
    - Np
    - Am
    Gd: []
    ```
    specifies that Gd will be the surrogate for no other elements (except of course itself), and Np and Am will both be treated as if they were U for the purposes of thermochemistry. There is one important keyword in this file: `Misc`, which is a catchall for all other elements that don't have defined surrogates. In the given `surrogates.yaml` these are all mapped to Pd, a noble metal, since they will benignly plate out and not influence the thermochemistry. Note that this list may include more elements (as keys) than `thermo_elements.csv`, but you must ensure that if no element with assigned to `Misc` is present, that your `thermo_elements.csv` list is inclusive of the elements of all isotopes in `isotopes.csv` otherwise surrogates are not able to be properly assigned for all elements.
## Full-Scale Parameters
The current parameters in the `simulation_parameters` (in terms of `isotopes.csv`, `species_elements.csv`, and `thermo_elements.csv`) are not representative of the actual model. They are a reduced set of parameters mainly for proofing the workflow of the simulation, and result in much faster running simulations for testing. The parameters designed for the final simulation are contained the the directory `source_term/data/simulation_parameters/full_scale` and must be replaced with those in the parent directory. This will result in a very long simulation simply due to the number of isotopes, and chemical/transport species being tracked. 

## Initial Isotopics
Currently the initial isotopics of the fuel salt are hardcoded based on a representative salt with impurities, but it may be changed to be more similar to the salt composition used for generating the cross sections by updating the `initial_isotopics` `defaultdict` in `input_generator.py`.

# Getting Started
## Add Common to PYTHONPATH
First add the directory `common` to your `PYTHONPATH`. For example, if this repository is cloned in `~/msr_source_term` add the following to your `.bashrc`
```
export PYTHONPATH=~/msr_source_term/common:$PYTHONPATH
```
## Install Python Requirements
Install the Python requirements (other than OpenMC) via
```
mamba install --file requirements.txt
```
or if using pip
```
pip install -r requirements.txt
```
## Edit Simulation Parameters
First, set the `blue_crab_executable` parameter in `source_term/data/simulation_parameters/parameters.yaml` to the path of your blue crab executable built with the proper experimental code and Thermochimica as described above. Also, please note the `moose-dev` environment (the name of the module you loaded using `ml`) you used to build your bluecrab executable with and set the `moose_dev_env` parameter to this in `source_term/data/simulation_parameters/parameters.yaml`.
## Run Prerequisite Simulations
To populate the data needed for the primary source term simulation, please follow the instructions in `keps_simulation` and `xs_gen` to run the simulations.