# Material properties fuel
rho = 3279. # density [kg / m^3]  (@1000K)
mu = 0.005926

################################################################################
# Species transport constants
################################################################################
Sc_t = 0.7
A_interface = 545.57 # Assumed from MCRE-ENG-PRSNT-0029Rev1A, (pi * (2 * r_pump_column)^2)
V_pump = 9930 # Obtained from paraview 
d = 0.10226 # m, pump hydraulic diameter

!include inputs/constants.i

################################################################################
# GEOMETRY
################################################################################
!include ../../inputs/th_mesh.i

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

# Global params
!include ../../inputs/th_global_params.i

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = INSFVRhieChowInterpolator
    u = superficial_vel_x
    v = superficial_vel_y
    w = superficial_vel_z
    a_u = a_u_var
    a_v = a_v_var
    a_w = a_w_var
    pressure = pressure
    block = 'reactor pipe pump'
  []
[]

# Variables
!include inputs/variables.i

# FV Kernels
!include inputs/fv_kernels.i

# User Objects
!include ../../inputs/th_user_objects.i

[ICs]
  [a_u_var]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = a_u_var
    from_variable = 'a_u_var'
  []
  [a_v_var]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = a_v_var
    from_variable = 'a_v_var'
  []
  [a_w_var]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = a_w_var
    from_variable = 'a_w_var'
  []
  [superficial_vel_x]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = superficial_vel_x
    from_variable = 'superficial_vel_x'
  []
  [superficial_vel_y]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = superficial_vel_y
    from_variable = 'superficial_vel_y'
  []
  [superficial_vel_z]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = superficial_vel_z
    from_variable = 'superficial_vel_z'
  []
  [pressure]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = pressure
    from_variable = 'pressure'
  []
  [TKE]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = TKE
    from_variable = 'TKE'
  []
  [TKED]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = TKED
    from_variable = 'TKED'
  []
  [T]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = T
    from_variable = 'T'
  []
  [mu_t]
    type = SolutionIC
    solution_uo = fluid_solution
    variable = mu_t
    from_variable = 'mu_t'
  []
[]

# Aux kernels
!include inputs/aux_kernels.i

# Aux variables
!include inputs/aux_variables.i
!include inputs/trivial_aux_variables.i
[AuxVariables]
  [a_u_var]
    type = MooseVariableFVReal
  []
  [a_v_var]
    type = MooseVariableFVReal
  []
  [a_w_var]
    type = MooseVariableFVReal
  []
  [superficial_vel_x]
    type = INSFVVelocityVariable
  []
  [superficial_vel_y]
    type = INSFVVelocityVariable
  []
  [superficial_vel_z]
    type = INSFVVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [TKE]
    type = INSFVEnergyVariable
  []
  [TKED]
    type = INSFVEnergyVariable
  []
  [T]
    type = INSFVEnergyVariable
  []
  [mu_t]
    type = MooseVariableFVReal
  []
  [fiss_rate_density]
    type = MooseVariableFVReal
  []
  [K] # Mass transfer coefficient
    type = MooseVariableFVReal
  []
[]

[FunctorMaterials]
  [turbulent_diffusion_coef]
    type = ADParsedFunctorMaterial
    expression = '(mu_t / ${rho}) / ${Sc_t}'
    functor_symbols = 'mu_t'
    functor_names = 'mu_t'
    property_name = 'turbulent_diffusion_coef'
  []
[]

# Postprocessors
!include inputs/postprocessors.i
!include inputs/trivial_postprocessors.i
!include inputs/element_postprocessors.i
!include inputs/trivial_element_postprocessors.i

[Postprocessors]
  [fiss_rate_density]
    type = ElementExtremeValue
    variable = fiss_rate_density
    execute_on = 'timestep_end'
  []
[]

[Controls]
  [solve_control]
    type = BoolFunctionControl
    function = solve_fn
    parameter = '*/*/solve'
    execute_on = 'timestep_begin'
  []
[]

# Functions
!include inputs/functions.i
[Functions]
  [solve_fn]
    type = ParsedFunction
    expression = 'if(t > ${fparse dt/th_num_steps}, 1, 0)'
  []
[]

# Executioner
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type  -pc_hypre_type'
  petsc_options_value = 'hypre  boomeramg'
  line_search = 'none'
  nl_abs_tol = '1e-20'
  nl_max_its = 30
  automatic_scaling = true
  compute_scaling_once = false
  num_steps = ${th_num_steps}
  end_time = ${dt}
[]

# Outputs
!include inputs/outputs.i