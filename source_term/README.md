# Overview
## The Main Simulation Script `run_dep.py`
The main simulation files for the depletion, thermochemistry, and species transport are `run_dep.i`, `sub_thermo.i`, and `species/sub_species.i` respectively. To run the source term simulation, you may run the `run_dep.py` script with the arguments `-p` for the number of MPI procs and the number of cores per proc. This will first generate the input files using the `input_generator.py` script, then it will use the bluecrab executable and moose dev environment specified in the simulation parameters to run the input. In order to ensure proper power normalization with this experimental depletion feature, the depletion can only be run for a constant power interval. So, for each constant power interval in the `power_history.csv` file, a depletion simulation is run with initial isotopics taken from the previous interval (if any, if not the initial isotopics). This python script is meant to automate this entire process, but if for some reason you would like to do this manually, you may run the script with the 'input only' option, specifying the index of the power history you'd like to generate inputs for (assuming you have run all of the previous indices and have the depletion outputs for those, which inform the initial isotopics of your current simulation). For example
```
python run_dep.py -io 0
```
will generate inputs for the initial power interval, and will use the initial isotopics, then
```
python run_dep.py -io 1
```
will generate inputs for the next interval assuming you have output isotopics from the simulation of the first interval, and so on.
## Data Visualization and Dose Conversion
Once you have run the main simulation, you may visualize the isotopics by running `plot_isotopes.py`. To convert the time series isotopics data to a final source term, there are two steps.
### Convert the Isotopics Activites
In doing this, there are two possible assumptions
1. **The off gas is mitigated**, i.e. all of the off gas that exits the primary loop (through the pump bowl) is held up until end of life in an off gas system, giving the fission products more time to decay to safer levels, only being released at the very end.
2. **The off gas is unmitigated**, i.e. all of the off gas that escapes the primary system is immediately vented to the environment with no additional time to decay. In this case, the nuclide densities are not just read from the isotopics (which include decay in the "off gas compartment") rather they are computed by integrating the removal rates
$$
A_i = -\lambda_i\sum_j\int_{t_j}^{t_{j+1}}\frac{dN_{i,core}}{dt}dt = \lambda_i\sum_j\int_{t_j}^{t_{j+1}}\lambda_{removal,i,j}N_{i,core}dt
$$
Now, since we only have the removal constants $\lambda_{removal,i,j}$ and the nuclide density in the core at each of a discrete number of timesteps, and it is either monotone increasing or decreasing within that time period, we have $\max\left(\frac{dN_{i,core}}{dt}\big|_{t=t_{j}}, \frac{dN_{i,core}}{dt}\big|_{t=t_{j+1}}\right) = \max(\lambda_{removal,i,j}N_{i,core,j}, \lambda_{removal,i,j+1}N_{i,core,j+1})$ and so
$$
A_i = \lambda_i\sum_j\max(\lambda_{removal,i,j}N_{i,core,j}, \lambda_{removal,i,j+1}N_{i,core,j+1})
$$


For each of these cases, run the corresponding script, either `python calculate_mitigated_activities.py` or `python calculate_unmitigated_activities.py`.

### Convert the Activities to Dose
`dose_conversion/inl_convert_activities.py`. This script uses a specific source term conversion spreadsheet distributed by INL for the LOTUS cell specifically, and will automatically add the isotopics data to the spreadsheet and read the results from it. Unfortunately, there is no way to programmatically update the calculations in the spreadsheet once the nuclides are added, so it is necessary to do this conversion in two stages
1. **Add the Nuclides** Run `python inl_convert_activities.py --add_nuclides`, then open the spreadsheel in Excel and save it once the values update.
2. **Extract the Results** Run `python inl_conver_activities.py --extract_results`, which will print out the results and output them to a file: `results.csv`.