/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                        GOX                                   */
/*                                                              */
/*           (c) 2015 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "PorousMediaBase.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<PorousMediaBase>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("temperature_var", "nonlinear variable names that holds temperature");
  params.addCoupledVar("pressure_var", "nonlinear variable names that holds pressure");
  params.addCoupledVar("gas_mixture", "nonlinear variable names: O2, N2, CO, CO2, He, Ar,  ...");

  params.addParam<std::vector<Real> >("molecular_weights", "in grams/mole");
  params.addParam<std::vector<Real> >("e_kb", "in Kelvin, is e/kb see L-J interaction potential");
  params.addParam<std::vector<Real> >("sigma", "in Angstroms, L-J  interaction potential");
  params.addParam<std::vector<Real> >("dipole_moment", "in Debye, Stockmayerinteraction potential");
  params.addParam<std::vector<Real> >("polarizability", "in Anstrom^3, ");

  /// input parameters for gas molecular diffusion coefficients in gas mixture
  params.addParam<Real>("diffusivity_of_gas", 1.0e-6, "diffusivity of molecules in air, m^2/s");
  params.addParam<Real>("diffusivity_scaling_factor", 1.0, "between 0 and 1");

  params.addParam<Real>("initial_porosity", 0.15, "Initial porosity of medium");
  params.addParam<Real>("initial_bulk_density", 2000.0, "bulk density of porous media in Kg/m^3");
  params.addParam<Real>("molecular_weight", 12.01e-3, "moelcular weight of graphite in kg/mol");

  //params.addParam<Real>("grain_density",       2160.0,  "grain density of graphite in Kg/m^3");

  MooseEnum graphite_type("IG-110 NBG-18", "IG-110");

  params.addParam<MooseEnum>("graphite_type", graphite_type, "graphite types included: IG-110 or NBG-18");
  params.addParam<Real>("reactive_surface_area", 0.3, "initial reactive surface area in m^2/m^3");
  params.addParam<Real>("rate_scaling_factor",   1.0,     "scaling factor for reaction rate");
  params.addParam<Real>("power_exponent", 1.0, "exponent for O2 concentration in rate law");
  params.addParam<Real>("gas_const", 8.31434, "Gas constant, in J/mol K");
  params.addParam<Real>("system_temperature", 298.15, "System temperature at simulation, K");
  params.addParam<Real>("system_pressure",101325.0 , "System pressure at simulation, Pascal");

  return params;
}

PorousMediaBase::PorousMediaBase(const InputParameters & parameters)
  :Material(parameters),

   _input_diffusivity(getParam<Real>("diffusivity_of_gas")),
   _diffusivity_scaling_factor(getParam<Real>("diffusivity_scaling_factor")),

   _input_porosity(getParam<Real>("initial_porosity")),
   _input_bulk_density(getParam<Real>("initial_bulk_density")),
   _molecular_weight(getParam<Real>("molecular_weight")),
//   _grain_density(getParam<Real>("grain_density")),

   _graphite_type(getParam<MooseEnum>("graphite_type")),
   _input_reactive_area(getParam<Real>("reactive_surface_area")),
   _rate_scaling_factor(getParam<Real>("rate_scaling_factor")),

   _power_exponent(getParam<Real>("power_exponent")),

   /// some intrinsic material properties

   _diffusivity_of_O2(declareProperty<Real>("diffusivity_of_O2")),
   _diffusivity_of_N2(declareProperty<Real>("diffusivity_of_N2")),
   _diffusivity_of_CO(declareProperty<Real>("diffusivity_of_CO")),
   _diffusivity_of_CO2(declareProperty<Real>("diffusivity_of_CO2")),
   _diffusivity_of_H2O(declareProperty<Real>("diffusivity_of_H2O")),
   _diffusivity_of_H2(declareProperty<Real>("diffusivity_of_H2")),
   _diffusivity_of_NH3(declareProperty<Real>("diffusivity_of_NH3")),
   _diffusivity_of_He(declareProperty<Real>("diffusivity_of_He")),
   _diffusivity_of_Ar(declareProperty<Real>("diffusivity_of_Ar")),

   _porosity(declareProperty<Real>("porosity")),
   _porosity_old(declarePropertyOld<Real>("porosity")),

   _bulk_density(declareProperty<Real>("bulk_density")),
   _bulk_density_old(declarePropertyOld<Real>("bulk_density")),

   _conversion_factor(declareProperty<Real>("conversion_factor")),

   _SA(declareProperty<Real>("reactive_surface_area")),
   _k_eff(declareProperty<Real>("effective_reaction_rate")),
   _CO_to_CO2_ratio(declareProperty<Real>("CO_to_CO2_ratio")),

   _has_temp(isCoupled("temperature_var")),
   _temperature(_has_temp ? coupledValue("temperature_var") : _zero),

   _has_pressure(isCoupled("pressure_var")),
   _pressure(_has_pressure ? coupledValue("pressure_var") : _zero),

   _gas_const(getParam<Real>("gas_const")),

   _system_temperature(getParam<Real>("system_temperature")),
   _system_pressure(getParam<Real>("system_pressure"))

//    mat(NULL),
//    rhs(NULL),
//    sol(NULL),
//    values(NULL),
//    rows(NULL),
//    cols(NULL)
{
  _grain_density = _input_bulk_density / (1.0 - _input_porosity);

  unsigned int n_mixture = coupledComponents("gas_mixture");
  nl_vars.resize(n_mixture);
  _vals.resize(n_mixture);
  _grad_vals.resize(n_mixture);

  /// once Hai tried to use petsc to solve a linear system

// #if LIBMESH_HAVE_PETSC
//     MatCreateBAIJ(PETSC_COMM_SELF, _vals.size(), _vals.size(), _vals.size(), PETSC_DETERMINE, PETSC_DETERMINE, _vals.size(),NULL, _vals.size(),NULL,&mat);
//     VecCreate(PETSC_COMM_SELF,&rhs);
//     VecCreate(PETSC_COMM_SELF,&sol);
//     VecSetSizes(rhs,n_mixture,PETSC_DETERMINE);
//     VecSetSizes(sol,n_mixture,PETSC_DETERMINE);
//     PetscCalloc1(n_mixture,&cols);
//     PetscCalloc1(n_mixture,&rows);
//     PetscCalloc1(n_mixture*n_mixture,&values);
// #endif

  for (unsigned int i = 0; i < n_mixture; ++i)
  {
    _vals[i] = &coupledValue("gas_mixture", i);
    _grad_vals[i] = &coupledGradient("gas_mixture", i);

    MooseVariable * var(getVar("gas_mixture", i));
    nl_vars[i] = var->name();
  }

  _molecular_weights = getParam<std::vector<Real>>("molecular_weights");
  _e_kb = getParam<std::vector<Real>>("e_kb");
  _sigma = getParam<std::vector<Real>>("sigma");
  _dipole_moment = getParam<std::vector<Real>>("dipole_moment");
  _polarizability = getParam<std::vector<Real>>("polarizability");
}

PorousMediaBase::~PorousMediaBase()
{
//  MatDestroy(&mat);
//  VecDestroy(&rhs);
//  VecDestroy(&sol);
//  PetscFree(values);
//  PetscFree(cols);
//  PetscFree(rows);
}

void
PorousMediaBase::initStatefulProperties(unsigned n_points)
{
  for (unsigned int qp = 0; qp < n_points; ++qp)
  {
    _porosity[qp] = _input_porosity;
    _porosity_old[qp] = _input_porosity;

    _bulk_density[qp] = _input_bulk_density;
    _bulk_density_old[qp] = _input_bulk_density;
  }
}

void
PorousMediaBase::computeProperties()
{
  for(unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    Real T;
    Real P;
    Real R = _gas_const;

    if (_has_temp)
      T = _temperature[qp];
    else
      T = _system_temperature;

    if (_has_pressure)
      P = _pressure[qp];
    else
      P = _system_pressure;

    /// Initialize the diffusivity of gas molecuales from input
    /// These diffusivities are acctually calculated later on using
    /// kinetic gas theory

    _diffusivity_of_O2[qp]  = _input_diffusivity;
    _diffusivity_of_N2[qp]  = _input_diffusivity;
    _diffusivity_of_CO[qp]  = _input_diffusivity;
    _diffusivity_of_CO2[qp] = _input_diffusivity;
    _diffusivity_of_H2O[qp] = _input_diffusivity;
    _diffusivity_of_H2[qp]  = _input_diffusivity;
    _diffusivity_of_NH3[qp] = _input_diffusivity;
    _diffusivity_of_He[qp]  = _input_diffusivity;
    _diffusivity_of_Ar[qp]  = _input_diffusivity;

    /// Initializing the binary diffusion coefficent calculation at given (P, T)

    unsigned int n_mixture = nl_vars.size();
    std::vector<std::vector<Real> > binary_diffusion(n_mixture);
    for (unsigned int i = 0; i < n_mixture; ++i)
      binary_diffusion[i].resize(n_mixture);

    /// Constants

    Real kbJ = 1.38065e-23; /// [J/K], Boltzmann constant
    Real kbE = 1.38065e-16; /// [Erg/K]
    Real NAvo = 6.0221409e23; /// Avogadro's constant

    /// TODO need a reference
    /// Calculating the  binary diffusion coefficent  at given (P, T)

    for (unsigned int i = 0; i < n_mixture; ++i)
    {
      for (unsigned int j = i; j < n_mixture; ++j)
      {
        Real tiny = 1.0e-10;

        Real s1 = _sigma[i];
        Real s2 = _sigma[j];

        Real e1 = _e_kb[i];
        Real e2 = _e_kb[j];

        Real u1 = _dipole_moment[i];
        Real u2 = _dipole_moment[j];

        bool nonpolar1(u1 <= tiny);
        bool nonpolar2(u2 <= tiny);

        Real s12;
        Real e12;
        Real dr12;

        if (nonpolar1 && nonpolar2)  /// both species are nonpolar
        {
          s12 = 0.5 * (s1+ s2);
          e12 = std::sqrt(e1 * e2);
          dr12 = 0.0;
        }
        else if (nonpolar1 && !nonpolar2)
        {
          Real alfa1 = _polarizability[i] / std::pow(s1,3.0);
          Real ur2 = u2 * 1.0e-18 / std::sqrt(e2 * kbE * std::pow(s2 * 1.0e-8, 3));
          Real ep = 1.0 + 0.25 * alfa1 * ur2 * std::sqrt(e2/e1);
          s12 = 0.5 * (s1 + s2) * std::pow(ep, -1.0/6.0);
          e12 = std::sqrt(e1 * e2) * std::pow(ep, 2.0);
          dr12 = 0.0;
        }
        else if (!nonpolar1 && nonpolar2)
        {
          Real alfa2 = _polarizability[j] / std::pow(s2,3.0);
          Real ur1 = u1 * 1.0e-18 / std::sqrt(e1 * kbE * std::pow(s1 * 1.0e-8, 3));
          Real ep = 1.0 + 0.25 * alfa2 * ur1 * std::sqrt(e1/e2);
          s12 = 0.5 * (s1 + s2) * std::pow(ep, -1.0/6.0);
          e12 = std::sqrt(e1 * e2) * std::pow(ep, 2.0);
          dr12 = 0.0;
        }
        else if (!nonpolar1 && !nonpolar2) //both species are polar
        {
          s12 = 0.5 * (s1 + s2);
          e12 = std::sqrt(e1 * e2);
          dr12 = 0.5 * (u1 * 1.0e-18 * u2 * 1.0e-18)/(e12 * kbE * std::pow(s12 * 1.0e-8, 3)); 
        }
        else
          mooseError("error in dipole moment combinations");

        /// Parameters for computing collision integral omiga(1,1) and omiga(2,2)

        Real a1 = 1.0548;
        Real a2 = 0.15504;
        Real a3 = 0.55909;
        Real a4 = 2.1705;
        Real a5 = 0.093193;
        Real a6 = 1.5;

        Real Tr = T / e12;
        Real f11 = 1.0 + (std::exp(a5/Tr) - std::exp(-a6/Tr)) * std::pow(dr12,2) / ( 2.0 + 2.5 * dr12);
        Real Om11 = (a1 * std::pow(Tr, -a2) + std::pow(Tr+a3, -a4)) * f11;

        Real M1 = _molecular_weights[i];
        Real M2 = _molecular_weights[j];
        Real m12 = M1 * M2 / (M1 + M2) / (1000.0 * NAvo) ;
        Real pi = 3.1415926;

        Real D12 = 3.0/16.0 * std::sqrt(2.0*pi*std::pow(kbJ*T,3.0)/m12) / (P*pi*std::pow(s12 *1.0e-10, 2.0)*Om11);

        binary_diffusion[i][j] = binary_diffusion[j][i] = D12;
      }
    }

    /// TODO need reference
    ///  effective diffusion coefficient based on composition of mixture and binary diffusion coefficient

    std::vector<Real> chi(n_mixture);  /// molar fraction of each species
    std::vector<Real> molar_conc(n_mixture);  /// molar concentrations of each species


    for (unsigned int i = 0; i < n_mixture; ++i)
    {
      if( (*_vals[i])[qp] >= 0.0)
        molar_conc[i] = (*_vals[i])[qp];
      else
        molar_conc[i] = 0.0; /// avoiding negative concentrations due to nuemrical oscillations
    }

    Real total_molar_conc(0.0);
    for (unsigned int i = 0; i < n_mixture; ++i)
      total_molar_conc += molar_conc[i];

    Real tiny = 1.0e-15;
    if (total_molar_conc <= tiny)
      mooseError("gas mixture near vacuum");

    /// molar fraction
    for (unsigned int i = 0; i < n_mixture; ++i)
      chi[i] = molar_conc[i]/total_molar_conc;

    /// salinity check chi(i) should adds up to 1
    Real sum_chi(0.0);
    for (unsigned int i = 0; i < n_mixture; ++i)
      sum_chi += chi[i];

    if (std::abs(sum_chi-1.0) >= tiny)
      mooseError("gas mixture composition error");

    Real O2_conc(0.0);

    for (unsigned int i = 0; i < n_mixture; ++i)
    {
      Real sum_chi_over_D(0.0);

      for (unsigned int j = 0; j < n_mixture; ++j)
      {
        if (j != i)
          sum_chi_over_D += chi[j] / binary_diffusion[i][j];
      }

      Real D_eff = (1.0 - chi[i]) / sum_chi_over_D;

      std::string name_mixture(nl_vars[i]);
      std::transform(name_mixture.begin(), name_mixture.end(), name_mixture.begin(), ::tolower);

      if (name_mixture == "o2")
      {
        /// Find the nonlinear variables that holds oxygen concentration
        O2_conc = molar_conc[i];

        _diffusivity_of_O2[qp] = D_eff;
      }
      else if (name_mixture == "n2")
      {
        _diffusivity_of_N2[qp] = D_eff;
      }
      else if (name_mixture == "co")
      {
        _diffusivity_of_CO[qp] = D_eff;
      }
      else if (name_mixture == "co2")
      {
        _diffusivity_of_CO2[qp] = D_eff;
      }
      else if (name_mixture == "h2o")
      {
        _diffusivity_of_H2O[qp] = D_eff;
      }
      else if (name_mixture == "h2")
      {
        _diffusivity_of_H2[qp] = D_eff;
      }
      else if (name_mixture == "nh3")
      {
        _diffusivity_of_NH3[qp] = D_eff;
      }
      else if (name_mixture == "he")
      {
        _diffusivity_of_He[qp] = D_eff;
      }
      else if (name_mixture == "ar")
      {
        _diffusivity_of_Ar[qp] = D_eff;
      }
      else
        mooseError("gas species not defined in material kernel");
    }

    /// TODO need reference
    /// Calculate reaction rate constant based on formula provided by Joshua

    Real Le = 2.01e-5;

    Real kA10 = 2.09e3;
    Real EA1 = 71000;
    Real kA1 = kA10 * std::exp(-EA1/(R*T));

    Real kA20 = 1.80e3;
    Real EA2 = 28000;
    Real kA2 = kA20 * std::exp(-EA2/(R*T));

    Real kS0 = 9.02e7;
    Real ES = 113000.0;
    Real kS = kS0 * std::exp(-ES/(R*T));

    Real kD10 = 8.21e8;
    Real ED1 = 200000.0;
    Real kD1 = kD10 * std::exp(-ED1/(R*T));

    Real kD20 = 3.31e8;
    Real ED2 = 148000.0;
    Real kD2 = kD20 * std::exp(-ED2/(R*T));

    Real kD30 = 4.46E16;
    Real ED3 = 289000;
    Real kD3 = kD30 * std::exp(-ED3/(R*T));

    /// Eq. (13) in J.J. Kane et al. Journal of Nuclear Materials 493 (2017) 343-367
    Real kTop = (kA1*kD1+(kA1+kA2)*(kS+kD1)) * (kD2 + kD3);
    Real kBot = (kA1+kA2) * (kS+kD1) * O2_conc + (kA1 * O2_conc + kS + kD1) * (kD2 + kD3);
    Real kconst = kTop * Le/kBot;

    /// TODO Eq. (19) in J.J. Kane et al. Carbon 59(2013) 49-64
    /// reaction product ratio between CO and CO2: CO/CO2
    _CO_to_CO2_ratio[qp] = kD3 * (kA1+kA2) * (kS+kD1) / (kA1*kD1*(kS+kD1)+kD2*(kA1+kA2)*(kS+kD1));

//    // compute reactive surface area as a function of porosity changes
//     Real init_porosity(_input_porosity);
//     Real last_porosity(_porosity[qp]);
//     Real  k_factor (1.5);
//     Real  con_tmp = k_factor * (1.0 - (1.0 - last_porosity)/(1.0 - init_porosity));
//     if (con_tmp > 1.0) con_tmp = 1.0;

    /// TODO need reference
    /// Conversion factor at the end of previous time step

    Real kalfa = 1.0 - _bulk_density_old[qp] / _input_bulk_density;
    if (kalfa > 1.0) kalfa = 1.0;
    if (kalfa < 0.0) kalfa = 0.0;

    Real init_area(_input_reactive_area);

    /// TODO polynomials unpublished; in preparation

    switch (_graphite_type)
    {
    case 0:  /// IG-110
      _SA[qp] = init_area * (- 2115.80 * std::pow(kalfa,6)
                             + 6464.50 * std::pow(kalfa,5)
                             - 7465.00 * std::pow(kalfa,4)
                             + 4094.40 * std::pow(kalfa,3)
                             - 1137.80 * std::pow(kalfa,2)
                             +  158.70 * kalfa
                             +    1.00 );
      break;

    case 1:  /// NBG-18
      _SA[qp] = init_area * (-  630.71 * std::pow(kalfa,6)
                             + 1816.10 * std::pow(kalfa,5)
                             - 1836.00 * std::pow(kalfa,4)
                             +  800.05 * std::pow(kalfa,3)
                             -  216.09 * std::pow(kalfa,2)
                             +   64.30 * kalfa
                             +    2.30 );
      break;

    default: mooseError("Unknown graphite type");
      break;
    }

    /// Kinteic reaction rate (should be always >= 0.0)
    Real kinetic_rate = _rate_scaling_factor * kconst * _SA[qp] * std::pow(O2_conc, _power_exponent);

    if (kinetic_rate < 0.0)
      kinetic_rate = 0.0;

    /// Simplified effectice rate provided by Joshua
    _k_eff[qp] = kinetic_rate;

    /// Current bulk density
    _bulk_density[qp] = _bulk_density_old[qp] - kinetic_rate * _dt;

    /// Fully_reacted case
    if (_bulk_density[qp] < 0.0)
      _bulk_density[qp] = 0.0;

    _conversion_factor[qp] = 1.0 - _bulk_density[qp] / _input_bulk_density;

    /// As bulk density changes due to reaction, so does porosity
    _porosity[qp] = _porosity_old[qp] + (_bulk_density_old[qp] - _bulk_density[qp]) / _grain_density;

    if (_porosity[qp] > 1.0)
      _porosity[qp] = 1.0;

    /// TODO need reference
    /// modifying the diffusion coefficients due to porosity/tortuosity of the porenetwork
    /// this could be as complicated as needed in future implemntations
    /// this was the old implemntation a while ago
    /// Real scaling_factor = (0.90585 + 0.05052 * _input_bulk_density/1000) * (kalfa - 1.0) + 1;

    Real Z(0.0);

    switch (_graphite_type)
    {
    case 0: /// IG-110
      Z = 5.95e-3;
      break;
    case 1: /// NBG-18
      Z = 1.13e-3;
      break;
    default: mooseError("Unknown graphite type");
      break;
    }

    Real scaling_factor =  (1.0 - Z) * kalfa + Z;

    _diffusivity_of_O2[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_O2[qp];
    _diffusivity_of_N2[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_N2[qp];
    _diffusivity_of_CO[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_CO[qp];
    _diffusivity_of_CO2[qp] = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_CO2[qp];
    _diffusivity_of_H2O[qp] = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_H2O[qp];
    _diffusivity_of_H2[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_H2[qp];
    _diffusivity_of_NH3[qp] = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_NH3[qp];
    _diffusivity_of_He[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_He[qp];
    _diffusivity_of_Ar[qp]  = _diffusivity_scaling_factor * scaling_factor * _diffusivity_of_Ar[qp];

// #if LIBMESH_HAVE_PETSC
//     for (PetscInt ii=0; ii<_vals.size(); ii++){
//       cols[ii] = ii;
//       rows[ii] = ii;
//       for(PetscInt jj=0; jj<_vals.size(); jj++)
//          if (ii==jj ) values[ii*_vals.size()+jj] = 1- chi[ii];
//          else  values[ii*_vals.size()+jj] = -chi[ii];
//     }
//     MatSetValues(mat,_vals.size(),rows, _vals.size(), cols, values, INSERT_VALUES);
//     MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
//     MatILUFactor(mat,NULL,NULL,NULL);
//     PetscScalar *array_rhs;
//     VecGetArray(rhs,&array_rhs);
//     for (PetscInt ii=0; ii<_vals.size(); ii++)
//       array_rhs[ii] = (*_grad_vals[ii])[qp](0);
//     VecRestoreArray(rhs,&array_rhs);
//     MatSolve(mat,rhs,sol);
//     PetscScalar *array_sol;
//     PetscScalar tmp;
//     VecGetArray(sol,&array_sol);
//     for (PetscInt ii=0; ii<_vals.size(); ii++)
//        tmp = array_sol[ii];
//     VecRestoreArray(sol,&array_sol);
// #endif
  }
}

