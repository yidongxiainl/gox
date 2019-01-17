/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GoxHeatConduction.h"

template<>
InputParameters validParams<GoxHeatConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Heat conduction");
  return params;
}

GoxHeatConduction::GoxHeatConduction(const InputParameters & parameters) :
    Diffusion(parameters),
    _kT(getMaterialProperty<Real>("thermal_conductivity"))
{
}

Real
GoxHeatConduction::computeQpResidual()
{
  return _kT[_qp] * Diffusion::computeQpResidual();
}

Real
GoxHeatConduction::computeQpJacobian()
{
  /// temperature-dependency of thermal conductivity not considered

  return _kT[_qp] * Diffusion::computeQpJacobian();
}
