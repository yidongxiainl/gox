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

#include "BulkDensityAux.h"

template<>
InputParameters validParams<BulkDensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
//   params.addParam<Real>("log_k", 0.0, "The equilibrium constant of the dissolution reaction");
//   params.addParam<Real>("reactive_surface_area", 0.1, "Specific reactive surface area in m^2/L solution");
//   params.addParam<Real>("ref_kconst", 6.456542e-8, "Kinetic rate constant in mol/m^2 s");
//   params.addParam<Real>("e_act", 2.91e4, "Activation energy, J/mol");
//   params.addParam<Real>("gas_const", 8.31434, "Gas constant, in J/mol K");
//   params.addParam<Real>("ref_temp", 298.15, "Reference temperature, K");
//   params.addParam<Real>("sys_temp", 298.15, "System temperature at simulation, K");

  params.addCoupledVar("oxygen_var","O2" "The name of nonlear variable that holds O2 concentration");
  
  params.addParam<std::vector<Real> >("sto_v", "The stochiometric coeff of reactant species");
  return params;
}

BulkDensityAux::BulkDensityAux(const InputParameters & parameters)
  :AuxKernel(parameters),
//    _log_k(getParam<Real>("log_k")),
//    _reactive_surface_area(getParam<Real>("reactive_surface_area")),
//    _ref_kconst(getParam<Real>("ref_kconst")),
//    _e_act(getParam<Real>("e_act")),
//    _gas_const(getParam<Real>("gas_const")),
//    _ref_temp(getParam<Real>("ref_temp")),
//    _sys_temp(getParam<Real>("sys_temp")),
//   _sto_v(getParam<std::vector<Real> >("sto_v")),
   
   _SA(getMaterialProperty<Real>("reactive_surface_area")),
   _k_eff(getMaterialProperty<Real>("effective_reaction_rate"))   

{
  unsigned int n = coupledComponents("oxygen_var");
  _vals.resize(n);
  for (unsigned int i = 0; i < _vals.size(); ++i)
    _vals[i] = &coupledValue("oxygen_var", i);
}

Real
BulkDensityAux::computeValue()
{
//   Real kconst = _ref_kconst * std::exp(-_e_act * (1.0 / _ref_temp - 1.0 / _sys_temp) / _gas_const);
//   Real omega = 1.0;

//   if (_vals.size())
//   {
//     for (unsigned int i = 0; i < _vals.size(); ++i)
//     {
//       if ((*_vals[i])[_qp] < 0.0)
//         omega *= std::pow(0.0, _sto_v[i]);
//       else
//         omega *= std::pow((*_vals[i])[_qp], _sto_v[i]);
//     }
//   }

//   Real saturation_SI = omega / std::pow(10.0, _log_k);
//   Real kinetic_rate = _reactive_surface_area * kconst * (1.0 - saturation_SI);

//   if (std::abs(kinetic_rate) <= 1.0e-12)
//     kinetic_rate =0.0;

  Real kinetic_rate = _k_eff[_qp] * _SA[_qp] *  (*_vals[0])[_qp];
  Real _bulk_density = _u_old[_qp] - kinetic_rate * _dt;

  // fully_reacted case
  if ( _bulk_density < 0.0)
     _bulk_density = 0.0;

  return  _bulk_density;
}
