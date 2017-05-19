/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GASMIXTUREDIFFUSION_H
#define GASMIXTUREDIFFUSION_H

#include "Diffusion.h"
#include "Material.h"

//Forward Declarations
class GasMixtureDiffusion;

template<>
InputParameters validParams<GasMixtureDiffusion>();


class GasMixtureDiffusion : public Diffusion
{
public:

  GasMixtureDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned _dim;
  const MaterialProperty<Real> & _diffusion_coefficient;
  const MaterialProperty<Real> * const _diffusion_coefficient_dT;
};

#endif //GASMIXTUREDIFFUSION_H
