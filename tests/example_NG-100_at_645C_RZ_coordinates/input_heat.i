### The PDFs to be solved are based on
### "J.J. Kane et al. Journal of Nuclear Materials 493 (2017) 343-367"
### The Eq. numbers below are referred from this article

################################################################################
[Problem]
  coord_type = RZ
[]
################################################################################
[Mesh]
  file = mesh.e

  ### mesh information:
  ### dim = 2
  ### nx = 100
  ### ny = 200
  ### xmax = 0.0127 # Test chamber radius
  ### ymax = 0.0254 # Length of test chamber
[]
################################################################################
[Variables]
  ### Molar concentration of oxygen gas [mol/m^3]
  [./O2]
    initial_condition = 0.0
  [../]

  ### Molar concentration of nitrogen gas [mol/m^3]
  [./N2]
    initial_condition = 13.1711
  [../]

  ### Molar concentration of carbon-monoxide gas [mol/m^3]
  [./CO]
    initial_condition = 0.0
  [../]

  ### Molar concentration of carbon-dioxide gas [mol/m^3]
  [./CO2]
    initial_condition = 0.0
  [../]

  ### Molar concentration of helium gas [mol/m^3]
  [./He]
    initial_condition = 1.0e-3
  [../]

  ### Temperature [K]
  [./T]
    initial_condition = 918.15
  [../]
[]
################################################################################
[AuxVariables]
  [./bulk_density]
    order = CONSTANT
    family = Monomial
  [../]

  [./porosity]
    order = CONSTANT
    family = Monomial
  [../]

  [./O2_diffusivity]
    order = CONSTANT
    family = Monomial
  [../]

  [./CO2_diffusivity]
    order = CONSTANT
    family = Monomial
  [../]

  [./CO_diffusivity]
    order = CONSTANT
    family = Monomial
  [../]

  [./k_eff]
    order = CONSTANT
    family = Monomial
  [../]

  [./conversion_factor]
    order = CONSTANT
    family = Monomial
  [../]

  [./co_co2_ratio]
    order = CONSTANT
    family = Monomial
  [../]

  [./Cp_C]
    order = CONSTANT
    family = Monomial
  [../]

  [./Cp_CO]
    order = CONSTANT
    family = Monomial
  [../]

  [./Cp_CO2]
    order = CONSTANT
    family = Monomial
  [../]

  [./Cp_O2]
    order = CONSTANT
    family = Monomial
  [../]

  [./deltaHrxn]
    order = CONSTANT
    family = Monomial
  [../]

  [./deltaHrxn_CO]
    order = CONSTANT
    family = Monomial
  [../]

  [./deltaHrxn_CO2]
    order = CONSTANT
    family = Monomial
  [../]

  [./heatSourceRate]
    order = CONSTANT
    family = Monomial
  [../]

  [./kT]
    order = CONSTANT
    family = Monomial
  [../]
[]
################################################################################
[Kernels]
  ### Overall reaction network C + (1-1/2x)O2 = xCO + (1-x)CO2
  ###   x: fractionation between CO and CO2 sto_coefficients for O2,CO and CO2
  ###   O2: -[1 - (1/2)x]
  ###   CO:  x
  ###   CO2  1-x

  #################################
  ### O2 equation; refer to Eq.(21)
  #################################

  [./O2_time]
    type = GasMixtureTimeDerivative
    variable = O2
  [../]

  [./O2_diffusion]
    type = GasMixtureDiffusion
    variable = O2
    diffusion_coefficient_name = diffusivity_of_O2
  [../]

  [./O2_source_sink]
    type = ReactionSourceSinkKernel
    variable = O2
    #sto_coeff = -0.875
  [../]

  ##################################
  ### CO2 equation; refer to Eq.(23)
  ##################################

  [./CO2_time]
    type = GasMixtureTimeDerivative
    variable = CO2
  [../]

  [./CO2_diffusion]
    type = GasMixtureDiffusion
    variable = CO2
    diffusion_coefficient_name = diffusivity_of_CO2
  [../]

  [./CO2_source_sink]
    type = ReactionSourceSinkKernel
    variable = CO2
    #sto_coeff = 0.75
  [../]


  #################################
  ### CO equation; refer to Eq.(22)
  #################################

  [./CO_time]
    type = GasMixtureTimeDerivative
    variable = CO
  [../]

  [./CO_diffusion]
    type = GasMixtureDiffusion
    variable = CO
    diffusion_coefficient_name = diffusivity_of_CO
  [../]

  [./CO_source_sink]
    type = ReactionSourceSinkKernel
    variable = CO
    #sto_coeff = 0.25
  [../]

  #################################
  ### N2 equation; refer to Eq.(24)
  #################################

  [./N2_time]
    type = GasMixtureTimeDerivative
    variable = N2
  [../]

  [./N2_diffusion]
    type = GasMixtureDiffusion
    variable = N2
    diffusion_coefficient_name = diffusivity_of_N2
  [../]

  #################################
  ### He equation; refer to Eq.(24)
  #################################

  [./He_time]
    type = GasMixtureTimeDerivative
    variable = He
  [../]

  [./He_diffusion]
    type = GasMixtureDiffusion
    variable = He
    diffusion_coefficient_name = diffusivity_of_He
  [../]

  #################################
  ### Heat transport-reaction equation; refer to Eq.(25)
  #################################

  [./heat_time]
    type = GoxHeatTimeDerivative
    variable = T
  [../]

  [./heat_conduction]
    type = GoxHeatConduction
    variable = T
  [../]

  [./heat_reaction]
    type = GoxHeatReaction
    variable = T
  [../]
[]
################################################################################
[AuxKernels]
  [./bulk_density]
    type = MaterialRealAux
    variable = bulk_density
    property = bulk_density
    execute_on = 'initial timestep_end'
  [../]

  [./porosity]
    type = MaterialRealAux
    variable = porosity
    property = porosity
    execute_on = 'initial timestep_end'
  [../]

  [./O2_diff]
    type = MaterialRealAux
    variable = O2_diffusivity
    property = diffusivity_of_O2
    execute_on = 'initial timestep_end'
  [../]

  [./CO2_diff]
    type = MaterialRealAux
    variable = CO2_diffusivity
    property = diffusivity_of_CO2
    execute_on = 'initial timestep_end'
  [../]

  [./CO_diff]
    type = MaterialRealAux
    variable = CO_diffusivity
    property = diffusivity_of_CO
    execute_on = 'initial timestep_end'
  [../]

  [./rate_const]
    type = MaterialRealAux
    variable = k_eff
    property = effective_reaction_rate
    execute_on = 'initial timestep_end'
  [../]

  [./conversion]
    type = MaterialRealAux
    variable = conversion_factor
    property = conversion_factor
    execute_on = 'initial timestep_end'
  [../]

  [./CO_to_CO2]
    type = MaterialRealAux
    variable = co_co2_ratio
    property = CO_to_CO2_ratio
    execute_on = 'initial timestep_end'
  [../]

  [./Cp_C]
    type = MaterialRealAux
    variable = Cp_C
    property = heat_capacity_of_C
    execute_on = 'initial timestep_end'
  [../]

  [./Cp_CO]
    type = MaterialRealAux
    variable = Cp_CO
    property = heat_capacity_of_CO
    execute_on = 'initial timestep_end'
  [../]

  [./Cp_CO2]
    type = MaterialRealAux
    variable = Cp_CO2
    property = heat_capacity_of_CO2
    execute_on = 'initial timestep_end'
  [../]

  [./Cp_O2]
    type = MaterialRealAux
    variable = Cp_O2
    property = heat_capacity_of_O2
    execute_on = 'initial timestep_end'
  [../]

  [./deltaHrxn]
    type = MaterialRealAux
    variable = deltaHrxn
    property = deltaHrxn
    execute_on = 'initial timestep_end'
  [../]

  [./deltaHrxn_CO]
    type = MaterialRealAux
    variable = deltaHrxn_CO
    property = deltaHrxn_CO
    execute_on = 'initial timestep_end'
  [../]

  [./deltaHrxn_CO2]
    type = MaterialRealAux
    variable = deltaHrxn_CO2
    property = deltaHrxn_CO2
    execute_on = 'initial timestep_end'
  [../]

  [./heatSourceRate]
    type = MaterialRealAux
    variable = heatSourceRate
    property = heat_source_rate
    execute_on = 'initial timestep_end'
  [../]

  [./kT]
    type = MaterialRealAux
    variable = kT
    property = thermal_conductivity
    execute_on = 'initial timestep_end'
  [../]
[]
################################################################################
[BCs]
  [./O2]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = O2
    final = 2.77
    rate = 2.0
  [../]

  [./CO2]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = CO2
    final = 0.0
    rate = 2.0
  [../]

  [./CO]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = CO
    final = 0.0
    rate = 2.0
  [../]

  [./N2]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = N2
    final = 10.4052
    rate = 2.0
  [../]

  [./He]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = He
    final = 1.0e-3
    rate = 2.0
  [../]

  [./T]
    type = ConvectiveFluxBC
    boundary = '2 4'
    variable = T
    final = 918.15
    rate = 2.0
  [../]

#  [./T]
#    type = DirichletBC
#    boundary = '2 4'
#    variable = T
#    value = 918.15
#  [../]
[]
################################################################################
[Materials]
  [./porous]
    type = PorousMediaBase
    block = 1

    gas_mixture  = 'O2 N2 CO CO2 He'
    molecular_weights = '31.9988 28.0134 28.0101 43.9987 4.0026'
    e_kb = '106.7 71.4 91.7 195.2 10.22'
    sigma = '3.467 3.798 3.69 3.941 2.551'
    dipole_moment = '0.0 0.0 0.122 0.0 0.0'
    polarizability = '1.562 1.71 1.9532 2.5070 0.2080'

    graphite_type = IG-110 # two types considered: IG-110 or NBG-18

    ### input parameters for IG-110 at 645 degree Celsus

    initial_porosity = 0.2114
    initial_bulk_density = 1779
    reactive_surface_area = 2.5e3
    rate_scaling_factor = 0.0105 # 0.0138
    power_exponent = 1.0

    ### input parameters for IG-110 at 744 degree Celsus

#    initial_porosity = 0.2114
#    initial_bulk_density = 1779
#    reactive_surface_area = 2.5e3
#    rate_scaling_factor = 0.016
#    power_exponent = 1.0

    ### input parameters for IG-110 at 564 degree Celsus

#    initial_porosity = 0.2114
#    initial_bulk_density = 1779
#    reactive_surface_area = 2.5e3
#    rate_scaling_factor = 0.018
#    power_exponent = 1.0

     ### input parameters for NBG-18 at 645 degree Celsus

#     initial_porosity = 0.1782
#     initial_bulk_density = 1851
#     diffusivity_scaling_factor = 2.19697
#     reactive_surface_area = 2.5e3
#     rate_scaling_factor = 0.00475
#     power_exponent = 1.0

     ### default temperature

#     system_temperature = 837.15 # 564 degree Celsus
     system_temperature = 918.15 # 645 degree Celsus
#     system_temperature = 1017.15 # 744 degree Celsus

    ### default pressure

    system_pressure = 101325.0

    temperature_var = 'T'

    ### input thermal conductivity
    thermal_conductivity = 2000.0
  [../]
[]
################################################################################
[Executioner]
  type = Transient
  #solve_type = 'PJFNK'
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_ls -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 201 cubic 0.7'

  dtmax = 600.0
  dtmin = 1.0e-10
  end_time = 216000

  ############################
  ### ONLY for regression test
  ############################
#  num_steps = 5

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1e-3
    percent_change = 0.12
  [../]

  l_max_its = 30
  l_tol = 1e-5

  nl_max_its = 10
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-12
[]
################################################################################
# [Postprocessors]
#   [./dofs]
#     type = NumDOFs
#   [../]
#   [./integral]
#     type = ElementL2Error
#     variable = temp
#     function = Combust
#   [../]
# []

[Postprocessors]
  [./prop_integral]
    type = ElementIntegralMaterialProperty
    mat_prop = bulk_density
    execute_on = 'initial timestep_end'
    outputs=csv
  [../]
[]
################################################################################
[Outputs]
  file_base = out-heat-IG-110_645C
  interval = 1
  exodus = true
  csv = true
  [./Console]
    type = Console
    output_linear = true
    output_nonlinear = true
    interval = 1
  [../]
[]
################################################################################
