 [Problem]
    coord_type = RSPHERICAL # XYZ RZ
 []

[Mesh]
# file = mesh.e

   type = GeneratedMesh
   dim = 1
   nx = 200
#    ny = 200
   xmax  =  0.05 # radius of the sphere (in meters)
#    ymax = 0.0254 # Test chamber radius
[]

[Variables]
  [./O2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0 
  [../]
  [./N2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 13.1711
  [../]
  [./CO]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./CO2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./He]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0e-3
  [../]
[]

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

[]

[Kernels]
#  Overall reaction network C + (1-1/2x)O2 = xCO + (1-x)CO2
#  , where x is the fractionation between CO and CO2
#  sto_coefficients for O2,CO and CO2
#  O2: -[1 - (1/2)x]
#  CO:  x
#  CO2  1-x
  [./O2_time]
    type = GasMixtureTimeDerivative
    variable = O2
  [../]
  [./O2_diffusiomn]
    type = GasMixtureDiffusion
    variable = O2
    diffusion_coefficient_name = diffusivity_of_O2
  [../]
  [./O2_source_sink]
    type = ReactionSourceSinkKernel
    variable = O2
#    sto_coeff = -0.875
  [../]

  [./CO2_time]
    type = GasMixtureTimeDerivative
    variable = CO2
  [../]
  [./CO2_diffusiomn]
    type = GasMixtureDiffusion
    variable = CO2
    diffusion_coefficient_name = diffusivity_of_CO2
  [../]
  [./CO2_source_sink]
    type = ReactionSourceSinkKernel
    variable = CO2
#    sto_coeff = 0.75
  [../]

  [./CO_time]
    type = GasMixtureTimeDerivative
    variable = CO
  [../]
  [./CO_diffusiomn]
    type = GasMixtureDiffusion
    variable = CO
    diffusion_coefficient_name = diffusivity_of_CO
  [../]
  [./CO_source_sink]
    type = ReactionSourceSinkKernel
    variable = CO
#    sto_coeff = 0.25
  [../]

  [./nitrogen_time]
    type = GasMixtureTimeDerivative
    variable = N2
  [../]
  [./nitrogen_diffusiomn]
    type = GasMixtureDiffusion
    variable = N2
    diffusion_coefficient_name = diffusivity_of_N2
  [../]
  [./He_time]
    type = GasMixtureTimeDerivative
    variable = He
  [../]
  [./He_diffusiomn]
    type = GasMixtureDiffusion
    variable = He
    diffusion_coefficient_name = diffusivity_of_He
  [../]

[]

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

[]


[BCs]


  [./O2]
#    type = DirichletBC
    type = ConvectiveFluxBC
    variable = O2
#    value = 2.77
     final = 2.77
     rate = 2.0
    boundary = right
  [../]
  [./CO2]
#    type = DirichletBC
    type = ConvectiveFluxBC
    variable = CO2
#    value = 0
     final = 0.0
     rate = 2.0
    boundary = right
  [../]

  [./CO]
#    type = DirichletBC
    type = ConvectiveFluxBC
    variable = CO
#    value = 0.0
     final = 0.0
     rate = 2.0
    boundary = right
  [../]

  [./N2]
#    type = DirichletBC
    type = ConvectiveFluxBC
    variable = N2
#    value = 10.4052
     final = 10.4052
     rate = 2.0
    boundary = right
  [../]
  [./He]
#    type = DirichletBC
    type = ConvectiveFluxBC
    variable = He
#    value = 1.0e-3
     final = 1.0e-3
     rate = 2.0
    boundary = right
  [../]
[]
[Materials]
  [./porous]
    type = PorousMediaBase
    block = 0

    gas_mixture  = 'O2 N2 CO CO2 He'
    molecular_weights =  '31.9988 28.0134 28.0101 43.9987 4.0026'
    e_kb =              '106.7 71.4 91.7 195.2 10.22'
    sigma =             '3.467 3.798 3.69 3.941 2.551'
    dipole_moment  =     '0.0 0.0 0.122 0.0 0.0'
    polarizability =     '1.562 1.71 1.9532 2.5070 0.2080'

    graphite_type = IG-110 # two types considered: IG-110 or NBG-18

# input parameters for IG-110 at 645 degree Celsus
    initial_porosity   = 0.2114
    initial_bulk_density    = 1779
    reactive_surface_area = 2.5e3
    rate_scaling_factor = 0.0105 # 0.0138
    power_exponent = 1.0

# input parameters for IG-110 at 744 degree Celsus
#      initial_porosity   = 0.2114
#      initial_bulk_density    = 1779
#      reactive_surface_area = 2.5e3
#      rate_scaling_factor = 0.016
#      power_exponent = 1.0

# input parameters for IG-110 at 564 degree Celsus
#      initial_porosity   = 0.2114
#      initial_bulk_density    = 1779
#      reactive_surface_area = 2.5e3
#      rate_scaling_factor = 0.018
#      power_exponent = 1.0


# input parameters for NBG-18 at 645 degree Celsus
#      initial_porosity   = 0.1782
#      initial_bulk_density    = 1851
##      diffusivity_scaling_factor = 2.19697
#      reactive_surface_area = 2.5e3
#      rate_scaling_factor = 0.00475
#      power_exponent = 1.0
    
#    sys_temp     = 837.15 # 564 degree Celsus
    sys_temp     = 918.15 # 645 degree Celsus
#    sys_temp     = 1017.15 # 744 degree Celsus

    sys_pressure = 101325.0
 
  [../]
[]


[Executioner]
  type       = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_ls -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 201 cubic 0.7'


  dtmax = 600.0
  dtmin = 1.0e-10
  end_time = 216000
  

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 100
    percent_change = 0.12
  [../]

  l_max_its  = 50 
  l_tol      = 1e-5
  nl_max_its = 10 
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-12
[]

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
  [../]
[]

[Outputs]
  file_base      = new-IG-110_645C
  interval       = 1
  exodus         = true
  csv            = true
  [./Console]
    type = Console
    linear_residuals    = 1
    nonlinear_residuals = 1
  [../]
[]
