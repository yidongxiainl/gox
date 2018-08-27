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
   //  params.addParam<Real>("molecular_weight", 21.01e-3, "stochiometric coefficients of minerals");
  //  Hai double checked this, should this be the molecular weight (in kg/mole) of carbo?n
    params.addParam<Real>("molecular_weight", 12.01e-3, "stochiometric coefficients of minerals");

  //Overall reaction network C + (1-1/2x)O2 = xCO + (1-x)CO2
  //, where x is the fractionation between CO and CO2
  // sto_coefficients for O2,CO and CO2
  // O2: -[1 - (1/2)x]
  // CO:  x
  // CO2  1-x
  // params.addRequiredParam<std::vector<Real> >("sto_coeff","stochiometric coefficients of reactant species");

  return params;
}

ReactionSourceSinkKernel::ReactionSourceSinkKernel(const InputParameters & parameters) :
    TimeDerivative(parameters),
    _porosity(getMaterialProperty<Real>("porosity")),
    _bulk_density(getMaterialProperty<Real>("bulk_density")),
    _bulk_density_old(getMaterialPropertyOld<Real>("bulk_density")),
    _CO_to_CO2_ratio(getMaterialProperty<Real>("CO_to_CO2_ratio")),
    _molecular_weight(getParam<Real>("molecular_weight"))

{
//   int n = coupledComponents("bulk_density_var");
//   _droh_dt.resize(n);

//   for (unsigned int i=0; i<_droh_dt.size(); ++i)
//     _droh_dt[i] = &coupledDot("bulk_density_var", i);
}

Real
ReactionSourceSinkKernel::computeQpResidual()
{
//   Real re = 0.0;
//   for (unsigned int i=0; i < _droh_dt.size(); ++i)
//     re += _sto_weight[i] * (*_droh_dt[i])[_qp] / (_molecular_weight * _porosity[_qp]) * _test[_i][_qp];

  Real re = 0.0;
  Real droh_dt = (_bulk_density[_qp] - _bulk_density_old[_qp]) / _dt;
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
  
  // re = sto_v * droh_dt / (_molecular_weight * _porosity[_qp]) * _test[_i][_qp];
  // Hai's new implemntation on Aug 23, 2018
  re = sto_v * droh_dt / _molecular_weight * _test[_i][_qp];
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
