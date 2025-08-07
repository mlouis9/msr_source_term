[Mesh]
  [gmg]
    type = GeneratedIDMeshGenerator
    dim = 2
  []
[]

!include inputs/thermo_constants.i

# Chemical composition
[ChemicalComposition]
  [thermo]
    thermofile = ${thermofile}
    output_phases = 'MSCL gas_ideal'
    elements = ${elements}
    output_element_phases = 'ALL'
    output_element_potentials = 'mu:Cl'
    tunit = K
    punit = atm
    munit = moles
    temperature = ${temperature}
    pressure = ${pressure}
    output_species_unit = moles
    is_fv = true
  []
[]


# Initial conditions
!include inputs/ics.i

# Used to postprocess input element values from transfers. Since elements are constant monomials on a 
# single-element mesh, ElementExtremeValue just gives the value of the variable
!include inputs/thermo_postprocessors.i

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = 'timestep_end'
  [element_output]
    type = CSV
    time_column = 'false'
    file_base = 'element_values'
  []
[]