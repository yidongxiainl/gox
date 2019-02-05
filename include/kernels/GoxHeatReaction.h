//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GOXHEATREACTION_H
#define GOXHEATREACTION_H

#include "Kernel.h"

// Forward Declaration
class GoxHeatReaction;

template <>
InputParameters validParams<GoxHeatReaction>();

class GoxHeatReaction : public Kernel
{
public:
  GoxHeatReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:

  const MaterialProperty<Real> & _heatSourceRate;
};
#endif // GOXHEATREACTION_H
