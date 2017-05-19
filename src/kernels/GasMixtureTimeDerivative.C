/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "GasMixtureTimeDerivative.h"


template<>
InputParameters validParams<GasMixtureTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}


GasMixtureTimeDerivative::GasMixtureTimeDerivative(const InputParameters & parameters) :
    TimeDerivative(parameters),
    _porosity(getMaterialProperty<Real>("porosity"))
{}

Real
GasMixtureTimeDerivative::computeQpResidual()
{
  return _porosity[_qp] * TimeDerivative::computeQpResidual();
}

Real
GasMixtureTimeDerivative::computeQpJacobian()
{
  return _porosity[_qp] * TimeDerivative::computeQpJacobian();
}
