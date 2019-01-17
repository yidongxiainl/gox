/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GOXHEATCONDUCTION_H
#define GOXHEATCONDUCTION_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class GoxHeatConduction;

template<>
InputParameters validParams<GoxHeatConduction>();


class GoxHeatConduction : public Diffusion
{
public:

  GoxHeatConduction(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:

  const MaterialProperty<Real> & _kT;
};

#endif //GOXHEATCONDUCTION_H
