/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "GoxHeatTimeDerivative.h"


template<>
InputParameters validParams<GoxHeatTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}


GoxHeatTimeDerivative::GoxHeatTimeDerivative(const InputParameters & parameters)
  :TimeDerivative(parameters),
   _bulk_density(getMaterialProperty<Real>("bulk_density")),
   _cp(getMaterialProperty<Real>("heat_capacity")),
   _dRhoCpT_dt(getMaterialProperty<Real>("heat_time_derivative"))
{
}

Real
GoxHeatTimeDerivative::computeQpResidual()
{
  /// This formulation only applies to first-order time derivative

  return _dRhoCpT_dt[_qp] * _test[_i][_qp];
}

Real
GoxHeatTimeDerivative::computeQpJacobian()
{
  /// The temperature-dependent change of (rho*Cp) is neglected

  return _bulk_density[_qp] * _cp[_qp] * TimeDerivative::computeQpJacobian();
}
