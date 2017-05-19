/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GASMIXTURETIMEDERIVATIVE_H
#define GASMIXTURETIMEDERIVATIVE_H

// MOOSE includes
#include "TimeDerivative.h"
#include "Material.h"

// Forward Declarations
class GasMixtureTimeDerivative;

template<>
InputParameters validParams<GasMixtureTimeDerivative>();

class GasMixtureTimeDerivative : public TimeDerivative
{
public:

  GasMixtureTimeDerivative(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:

  const MaterialProperty<Real> & _porosity;
  
};

#endif //GASMIXTURETIMEDERIVATIVE_H
