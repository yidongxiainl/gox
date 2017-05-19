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

  
  //input parameters for gas molecular diffusion coefficients in gas mixture   
 params.addParam<Real>("diffusivity_of_gas",  1.0e-6, "diffusivity of molecules in air, m^2/s");
  params.addParam<Real>("diffusivity_scaling_factor",  1.0, " between 0 and 1");
  
  params.addParam<Real>("initial_porosity",         0.15,  "Initial porosity of medium");
  params.addParam<Real>("initial_bulk_density",       2000.0,  "bulk density of porous media in Kg/m^3");
  params.addParam<Real>("molecular_weight",      21.01e-3,  "moelcular weight of graphite in kg/mol");
//  params.addParam<Real>("grain_density",       2160.0,  "grain density of graphite in Kg/m^3");


  params.addParam<Real>("reactive_surface_area",   0.3,     "initial reactive surface area in m^2/m^3");
  params.addParam<Real>("rate_scaling_factor",   1.0,     "scaling factor for reaction rate");
  
  params.addParam<Real>("power_exponent", 1.0, "exponent for O2 concentration in rate law");


  params.addParam<Real>("gas_const", 8.31434, "Gas constant, in J/mol K");

  params.addParam<Real>("sys_temp", 298.15, "System temperature at simulation, K");
  params.addParam<Real>("sys_pressure",101325.0 , "System pressure at simulation, Pascal");
  
  return params;

}

PorousMediaBase::PorousMediaBase(const InputParameters & parameters)
  :Material(parameters),

   _input_default_diffusivity(getParam<Real>("diffusivity_of_gas")),
   _diffusivity_scaling_factor(getParam<Real>("diffusivity_scaling_factor")),

   _input_porosity(getParam<Real>("initial_porosity")),
   _input_bulk_density(getParam<Real>("initial_bulk_density")),
   _molecular_weight(getParam<Real>("molecular_weight")),
//   _grain_density(getParam<Real>("grain_density")),

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

   _sys_temp(getParam<Real>("sys_temp")),
   _sys_pressure(getParam<Real>("sys_pressure")),
   mat(NULL),
   rhs(NULL),
   sol(NULL),
   values(NULL),
   rows(NULL),
   cols(NULL)
{
  _grain_density = _input_bulk_density / (1.0 - _input_porosity);
  
  unsigned int n = coupledComponents("gas_mixture");
  nl_vars.resize(n);
  _vals.resize(n);
  _grad_vals.resize(n);
 
#if LIBMESH_HAVE_PETSC
    MatCreateBAIJ(PETSC_COMM_SELF, _vals.size(), _vals.size(), _vals.size(), PETSC_DETERMINE, PETSC_DETERMINE, _vals.size(),NULL, _vals.size(),NULL,&mat);
    VecCreate(PETSC_COMM_SELF,&rhs);
    VecCreate(PETSC_COMM_SELF,&sol);
    VecSetSizes(rhs,n,PETSC_DETERMINE);
    VecSetSizes(sol,n,PETSC_DETERMINE);
    PetscCalloc1(n,&cols);
    PetscCalloc1(n,&rows);
    PetscCalloc1(n*n,&values);
#endif 
  
  for (unsigned int i = 0; i < _vals.size(); ++i)
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
 MatDestroy(&mat);
 VecDestroy(&rhs);
 VecDestroy(&sol);
 PetscFree(values);
 PetscFree(cols);
 PetscFree(rows);
}

void
PorousMediaBase::initStatefulProperties(unsigned n_points)
{
  for (unsigned int qp = 0; qp < n_points; ++qp)
  {
    _porosity[qp] = _input_porosity;
    _bulk_density[qp] = _input_bulk_density;

    _porosity_old[qp] = _input_porosity;
    _bulk_density_old[qp] = _input_bulk_density;
  }
}

void
PorousMediaBase::computeProperties()
{
  for(unsigned int qp=0; qp<_qrule->n_points(); ++qp)
  {
    Real T;
    Real P;
    Real R = _gas_const;
    
    if (_has_temp)
      T = _temperature[qp];
    else
      T = _sys_temp;
    
    if (_has_pressure)
      P = _pressure[qp];
    else
      P = _sys_pressure;

    // Initialize the diffusivity of gas molecuales from input
    // These diffusivities are acctually calculated later on using
    // kinetic gas theory
    _diffusivity_of_O2[qp]  = _input_default_diffusivity;
    _diffusivity_of_N2[qp]  = _input_default_diffusivity;
    _diffusivity_of_CO[qp]  = _input_default_diffusivity;
    _diffusivity_of_CO2[qp] = _input_default_diffusivity;
    _diffusivity_of_H2O[qp] = _input_default_diffusivity;
    _diffusivity_of_H2[qp]  = _input_default_diffusivity;
    _diffusivity_of_NH3[qp] = _input_default_diffusivity;
    _diffusivity_of_He[qp]  = _input_default_diffusivity;
    _diffusivity_of_Ar[qp]  = _input_default_diffusivity;

    // Initializing the binary diffusion coefficent calculation at given (P, T)
    std::vector<std::vector<Real> > binary_diffusion(nl_vars.size());
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
      binary_diffusion[i].resize(nl_vars.size());

    //Constants
    Real kbJ=1.38065e-23; //in J/K, Boltzmann constant
    Real kbE=1.38065e-16; //in Erg/K
    Real NAvo=6.0221409e23; // Avogadro's constant

    // Calculating the  binary diffusion coefficent  at given (P, T)
//     _console << "*******************************************************************************************"<<std::endl;
//     _console << "                    ";
//     for (unsigned int i = 0; i < nl_vars.size(); ++i)
//          _console << nl_vars[i] + "              ";
//     _console << std::endl;
    
        
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
    {

      for (unsigned int j = i; j < nl_vars.size(); ++j)
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
        
        if (nonpolar1 && nonpolar2)  //both species are nonpolar
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
        else if (!nonpolar1 && !nonpolar2) //both species are po;ar
        {
          s12 = 0.5 * (s1 + s2);
          e12 = std::sqrt(e1 * e2);
          dr12 = 0.5 * (u1 * 1.0e-18 * u2 * 1.0e-18)/(e12 * kbE * std::pow(s12 * 1.0e-8, 3));
          
        }
        else
          mooseError("error in dipole moment combinations");

        //Parameters for computing collision integral omiga(1,1) and omiga(2,2)
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
        
        Real D12 = 3.0 / 16.0 *std::sqrt(2.0 * pi * std::pow(kbJ * T,3.0) / m12)
                       / (P * pi * std::pow(s12 *1.0e-10, 2.0) * Om11);

        binary_diffusion[i][j] = binary_diffusion[j][i] = D12;

      
      }
//       _console << std::endl;
//       _console << nl_vars[i] + "               ";
//       for (unsigned int j = 0; j < nl_vars.size(); ++j)
//         _console << binary_diffusion[i][j] << "    ";
//       _console << std::endl;
    }
    
    //  effective diffusion coefficient based on composition of mixture and binary diffusion coefficient
    std::vector<Real> X(nl_vars.size());  //mole fraction of each species
    std::vector<Real> molar_conc(nl_vars.size());  //molar concentrations of each species

    Real tiny = 1.0e-15;

    for (unsigned int i = 0; i < nl_vars.size(); ++i)
    {
      if( (*_vals[i])[qp] >= 0.0)
        molar_conc[i] = (*_vals[i])[qp];
      else
        molar_conc[i] = 0.0; //avoiding negative concentrations due to nuemrical oscillations
    }
    
    
    Real total_molar_conc(0.0);
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
      total_molar_conc += molar_conc[i];

    if (total_molar_conc <= tiny)
      mooseError("gas mixture near vacuum");

    // mole fraction
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
      X[i] = molar_conc[i]/total_molar_conc;

    // salinity check X(i) should adds up to 1
    Real sumX(0.0);
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
      sumX += X[i];

    if (std::abs(sumX - 1.0) >= tiny)
      mooseError("gas mixture composition error");

    for (unsigned int i = 0; i < nl_vars.size(); ++i)
    {

      Real sum_X_over_D(0.0);
      for (unsigned int j = 0; j < nl_vars.size(); ++j)
      {
        if (j != i)
          sum_X_over_D += X[j] / binary_diffusion[i][j];
      }

      Real D_eff = (1.0 - X[i]) / sum_X_over_D;

      std::string n(nl_vars[i]);
      std::transform( n.begin(), n.end(), n.begin(), ::tolower);
      if (n == "o2")
      {
        _diffusivity_of_O2[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      else if (n == "n2")
      {
        _diffusivity_of_N2[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "co")
      {
        _diffusivity_of_CO[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "co2")
      {
        _diffusivity_of_CO2[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      else if (n == "h2o")
      {
        _diffusivity_of_H2O[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "h2")
      {
        _diffusivity_of_H2[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "nh3")
      {
        _diffusivity_of_NH3[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "he")
      {
        _diffusivity_of_He[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else if (n == "ar")
      {
        _diffusivity_of_Ar[qp] = D_eff;
//        _console << n << " =  " << D_eff << std::endl;
      }
      
      else
        mooseError("gas species not defined in material kernel");
    }
    

    //find the nonlinear variables that holds oxygen concentration
    Real O2_conc(0.0);
    for (unsigned int i = 0; i < nl_vars.size(); ++i)
    {
      std::string n(nl_vars[i]);
      std::transform( n.begin(), n.end(), n.begin(), ::tolower);
      if (n == "o2")
        O2_conc = molar_conc[i];
    }
    //_console << O2_conc << "  "  << R <<"  "<< T <<std::endl;
    

    //Calculate reaction rate constant based on Josh's formular
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

    Real kTop = (kA1*kD1+(kA1+kA2)*(kS+kD1)) * (kD2 + kD3);
    Real kBot = (kA1+kA2) * (kS+kD1) * O2_conc + (kA1 * O2_conc + kS + kD1) * (kD2 + kD3);

    Real kconst = kTop * Le/kBot;

    // reaction product ratio between CO and CO2: CO/CO2
    Real Ra = kD3 * (kA1 + kA2) * (kS+kD1)/(kA2*kD2*(kS+kD1)+kD2*(kA1+kA2)*(kS+kD1));
    
    
    _CO_to_CO2_ratio[qp] = Ra;

//    _console << " effective reaction rate "<< kconst << "   "<< Ra << std::endl;

    // compute reactive surface area as a function of porosity changes
    Real init_area(_input_reactive_area);
    
    Real init_porosity(_input_porosity);
    Real last_porosity(_porosity[qp]);
    Real  con_tmp = 1.0 - (1.0 - last_porosity)/(1.0 - init_porosity);
    
//    _SA[qp]    = _input_reactive_area;
    _SA[qp]    = init_area * (  27.72 * std::pow(con_tmp,3)
                              - 77.04 * std::pow(con_tmp,2)
                              + 48.67 * con_tmp
                              +  1.0) ;
    
// kinteic reaction rate
    Real kinetic_rate = _rate_scaling_factor * kconst * _SA[qp] *  std::pow(O2_conc, _power_exponent); //Notice: this rate is always >= 0.0
    
    if (kinetic_rate < 0.0)
      kinetic_rate = 0.0;

    _k_eff[qp] = kinetic_rate; //this is based on Josh's simple effectice rate concept
    
    _bulk_density[qp] = _bulk_density_old[qp] - kinetic_rate * _dt;
    // fully_reacted case
    if ( _bulk_density[qp] < 0.0)
      _bulk_density[qp] = 0.0;

    _conversion_factor[qp] = 1.0 - _bulk_density[qp] / _input_bulk_density;
    
    // As bulk density changes due to reaction, so does porosity
    Real vf = ( _bulk_density_old[qp] -  _bulk_density[qp]) / _grain_density;
    _porosity[qp] = _porosity_old[qp] + vf;

    if (_porosity[qp] > 1.0)
      _porosity[qp] = 1.0;

    // modifying the diffusion coefficients due to porosity/tortuosity of the porenetwork
    // this could be as complicated as needed in future implemntations
    _diffusivity_of_O2[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_O2[qp];
    _diffusivity_of_N2[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_N2[qp];
    _diffusivity_of_CO[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_CO[qp];
    _diffusivity_of_CO2[qp] = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_CO2[qp];
    _diffusivity_of_H2O[qp] = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_H2O[qp];
    _diffusivity_of_H2[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_H2[qp];
    _diffusivity_of_NH3[qp] = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_NH3[qp];
    _diffusivity_of_He[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_He[qp];
    _diffusivity_of_Ar[qp]  = _diffusivity_scaling_factor * (1.21585 * _porosity[qp] - 0.215854) * _diffusivity_of_Ar[qp];

#if LIBMESH_HAVE_PETSC
    for (PetscInt ii=0; ii<_vals.size(); ii++){
      cols[ii] = ii;
      rows[ii] = ii;
      for(PetscInt jj=0; jj<_vals.size(); jj++)
         if (ii==jj ) values[ii*_vals.size()+jj] = 1- X[ii];
         else  values[ii*_vals.size()+jj] = -X[ii];
    }
    MatSetValues(mat,_vals.size(),rows, _vals.size(), cols, values, INSERT_VALUES);  
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    MatILUFactor(mat,NULL,NULL,NULL);
    PetscScalar *array_rhs;
    VecGetArray(rhs,&array_rhs);
    for (PetscInt ii=0; ii<_vals.size(); ii++)
      array_rhs[ii] = (*_grad_vals[ii])[qp](0);
    VecRestoreArray(rhs,&array_rhs);
    MatSolve(mat,rhs,sol); 
    PetscScalar *array_sol;
    PetscScalar tmp;
    VecGetArray(sol,&array_sol);
    for (PetscInt ii=0; ii<_vals.size(); ii++)
       tmp = array_sol[ii];
    VecRestoreArray(sol,&array_sol);  
#endif   
  }
}

