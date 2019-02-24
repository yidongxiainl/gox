/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GOXHEATTIMEDERIVATIVE_H
#define GOXHEATTIMEDERIVATIVE_H

// MOOSE includes
#include "TimeDerivative.h"
#include "Material.h"

// Forward Declarations
class GoxHeatTimeDerivative;

template<>
InputParameters validParams<GoxHeatTimeDerivative>();

class GoxHeatTimeDerivative : public TimeDerivative
{
public:

  GoxHeatTimeDerivative(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:

  const MaterialProperty<Real> & _bulk_density;
  const MaterialProperty<Real> & _cp_C;
  const MaterialProperty<Real> & _dRhoCpT_dt;
};

#endif //GOXHEATTIMEDERIVATIVE_H
