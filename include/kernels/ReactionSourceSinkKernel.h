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

#ifndef REACTIONSOURCESINKKERNEL
#define REACTIONSOURCESINKKERNEL

#include "TimeDerivative.h"

// Forward Declaration
class ReactionSourceSinkKernel;


template<>
InputParameters validParams<ReactionSourceSinkKernel>();


class ReactionSourceSinkKernel : public TimeDerivative
{
public:

  ReactionSourceSinkKernel(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();


  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Material property of porosity
  const MaterialProperty<Real> & _drho_dt;
  const MaterialProperty<Real> & _CO_to_CO2_ratio;

  Real _molecular_weight;

  /// Coupled time derivatives of mineral concentrations (stored and computed as Aux variables).
//  std::vector<const VariableValue *> _droh_dt;
};

#endif // REACTIONSOURCESINKKERNEL
