//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GoxHeatReaction.h"

template <>
InputParameters
validParams<GoxHeatReaction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Gox thermal reaction");
  return params;
}

GoxHeatReaction::GoxHeatReaction(const InputParameters & parameters) : Kernel(parameters) {}

Real
GoxHeatReaction::computeQpResidual()
{
  //return _test[_i][_qp] * _u[_qp];
  return 0.0;
}

Real
GoxHeatReaction::computeQpJacobian()
{
  //return _test[_i][_qp] * _phi[_j][_qp];
  return 0.0;
}
