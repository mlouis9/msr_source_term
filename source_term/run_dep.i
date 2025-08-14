# ---------
# Constants
# ---------
V =  331000                                       # cm3, from MCRE-ENG-PRSNT-0029 Rev1A
at_per_bcm_to_moles_per_cm3 = 1.6605390561442868  # 10^24 ( at/b*cm )/( at/cm^3 ) * 1/(6.0221408 x 10^23)( moles/at ) = 10/6.0221408
stability_eps = 1e-03                             # Excess chlorine (moles) to ensure Thermochimica GEM solve is stable
ar_moles = 2.021
max_removal_rate = 1

[Mesh]
  [gmg]
    type = GeneratedIDMeshGenerator
    dim = 2
    depletion_ids = 1
    material_ids = 1
  []
[]

[Problem]
  type = FEProblem
  solve = false
  verbose_multiapps = true
  verbose_setup = true
  kernel_coverage_check = false
[]

# Thermochimica AuxVariables
!include inputs/element_aux_variables.i
# Species Transport AuxVariables
!include inputs/th_aux_variables.i 
[AuxVariables]
  [tfuel]
    order = CONSTANT
    family = MONOMIAL
  []
  [flux_g0]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1.0
  []
  [depletion_s]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
  [fiss_rate_density]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
[]

!include inputs/element_aux_kernels.i
# Species transport AuxKernels
!include inputs/th_aux_kernels.i
[AuxKernels]
  [SetDepletion_s]
    type = FunctionAux
    function = SetDepletion_s
    variable = depletion_s
    execute_on = 'timestep_begin'
  []
  [set_fission_rate_density]
    type = FunctionAux
    variable = fiss_rate_density
    function = fiss_rate_density
    execute_on = 'timestep_begin'
  []
[]

# Executioner
!include inputs/dep_executioner.i

# Postprocessors
!include inputs/isotope_postprocessors.i
!include inputs/sumset_postprocessors.i
!include inputs/total_removal_postprocessors.i
# Species transport Postprocessors
!include inputs/th_postprocessors.i
# Vectorpostprocessors/bateman
!include inputs/vector_postprocessors.i

[Functions]
  [SetDepletion_s]
    type = PiecewiseConstant
    data_file = timesteps.csv
    format = 'columns'
    x_title = 'Timesteps'
    y_title = 'Depletion Values'
    direction = 'right_inclusive'
  []
  [PowerModulate]
    type = PiecewiseConstant
    data_file = 'data/simulation_parameters/power_history.csv'
    format = 'columns'
    x_title = 'Timesteps'
    y_title = 'Power Factor'
    direction = 'right_inclusive'
  []
  [fiss_rate_density]
    type = PiecewiseConstant
    data_file = 'data/simulation_parameters/fiss_rate_densities.csv'
    format = 'columns'
    x_title = 'Timesteps'
    y_title = 'Fission Rate Density'
    direction = 'right_inclusive'
  []
[]

# Total
!include inputs/total_removal_functions.i

# Outputs
!include inputs/outputs.i

# ----------------------
# Multiphysics Coupling
# ----------------------
# Controls
# !include inputs/dep_controls.i

# MultiApps
!include inputs/multiapps.i

# Thermochimica transfers
!include inputs/thermo_transfers.i
# Species transport transfers
!include inputs/th_transfers.i
!include inputs/th_thermo_removal_transfers.i
[Transfers]
  [fiss_rate_density]
    type = MultiAppGeometricInterpolationTransfer
    to_multi_app = species
    source_variable = fiss_rate_density
    variable = fiss_rate_density
  []
[]