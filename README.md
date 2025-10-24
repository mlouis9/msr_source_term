# Introduction

# Directory description
1. [keps_simulation](keps_simulation)
2. [source_term](source_term)
3. [xs_gen](xs_gen)

# Requirements
- **Python libraries**:
  - Pandas
- **OpenMC**: A working installation of OpenMC (with the nuclear data libraries) is required for generating the cross sections in `xs_gen`, guidance is given [here](https://docs.openmc.org/en/stable/quickinstall.html).
- The latest version of **MSTDB-TC** (Molten Salt Thermal Properties Database - Thermo Chemical). Can be obtained [https://mstdb.ornl.gov/](here). Then place the FactSage `.dat` file for the Chloride database (e.g. `MSTDB-TC_V4.0_Chlorides_No_Func.dat`) in `source_term/data/thermo`.

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