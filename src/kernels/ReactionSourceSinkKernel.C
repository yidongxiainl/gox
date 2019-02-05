/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                        GOX                                   */
/*                                                              */
/*           (c) 2016 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ReactionSourceSinkKernel.h"
#include "Material.h"
#include "MooseVariable.h"

template<>
InputParameters validParams<ReactionSourceSinkKernel>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<Real>("molecular_weight", 12.01e-3, "stochiometric coefficients of minerals");

  return params;
}

ReactionSourceSinkKernel::ReactionSourceSinkKernel(const InputParameters & parameters)
  :TimeDerivative(parameters),

   _drho_dt(getMaterialProperty<Real>("bulk_density_time_derivative")),
   _CO_to_CO2_ratio(getMaterialProperty<Real>("CO_to_CO2_ratio")),
   _molecular_weight(getParam<Real>("molecular_weight"))
{
}

Real
ReactionSourceSinkKernel::computeQpResidual()
{
  Real re = 0.0;
  std::string n(_var.name());

  std::transform( n.begin(), n.end(), n.begin(), ::tolower);

  Real Ra(_CO_to_CO2_ratio[_qp]);

  Real sto_v;

  Real X = Ra / (1.0 + Ra);

  if (n == "o2")
    sto_v = -1.0 * (1.0 - 1.0/2.0 * X);
  else if (n == "co")
    sto_v = X;
  else if (n == "co2")
    sto_v = 1.0 - X;
  else
    mooseError("reaction source/sink kernel only acts on O2, CO, and CO2");

  re = sto_v * _drho_dt[_qp] / _molecular_weight * _test[_i][_qp];
  return re;
}

Real
ReactionSourceSinkKernel::computeQpJacobian()
{
  return 0.0;
}

Real ReactionSourceSinkKernel::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0.0;
}
