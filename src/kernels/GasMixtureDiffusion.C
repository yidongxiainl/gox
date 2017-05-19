/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GasMixtureDiffusion.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<GasMixtureDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("Compute diffusion of gas molecules");
  params.addParam<MaterialPropertyName>("diffusion_coefficient_name",
                                        "diffusivity_of_O2",
                                        "Property name of gas diffusivity (Default: diffusivity_of_O2)");
  params.addParam<MaterialPropertyName>("diffusion_coefficient_dT_name",
                                        "diffusivity_of_O2_dT",
                                        "Property name of the derivative of the diffusivity with respect "
                                        "to the variable (Default: diffusivity_of_O2_dT)");
  return params;
}

GasMixtureDiffusion::GasMixtureDiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _dim(_subproblem.mesh().dimension()),
    _diffusion_coefficient(getMaterialProperty<Real>("diffusion_coefficient_name")),
    _diffusion_coefficient_dT(hasMaterialProperty<Real>("diffusion_coefficient_dT_name") ?
                              &getMaterialProperty<Real>("diffusion_coefficient_dT_name") : NULL)
{
}

Real
GasMixtureDiffusion::computeQpResidual()
{
  return _diffusion_coefficient[_qp] * Diffusion::computeQpResidual();
}

Real
GasMixtureDiffusion::computeQpJacobian()
{
  Real jac = _diffusion_coefficient[_qp] * Diffusion::computeQpJacobian();
  if (_diffusion_coefficient_dT)
    jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  return jac;
}
