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


GoxHeatTimeDerivative::GoxHeatTimeDerivative(const InputParameters & parameters) :
    TimeDerivative(parameters),
    _bulk_density(getMaterialProperty<Real>("bulk_density")),
    _cp(getMaterialProperty<Real>("heat_capacity"))
{}

Real
GoxHeatTimeDerivative::computeQpResidual()
{
  return _bulk_density[_qp] * _cp[_qp] * TimeDerivative::computeQpResidual();
}

Real
GoxHeatTimeDerivative::computeQpJacobian()
{
  /// the temperature-dependent change of (rho*Cp) is neglected

  return _bulk_density[_qp] * _cp[_qp] * TimeDerivative::computeQpJacobian();
}
