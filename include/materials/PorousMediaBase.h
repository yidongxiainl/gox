/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                        GOX                                   */
/*                                                              */
/*           (c) 2015 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef POROUSMEDIABASE_H
#define POROUSMEDIABASE_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

//Forward Declarations
class PorousMediaBase;

template<>
InputParameters validParams<PorousMediaBase>();

class PorousMediaBase : public Material
{
public:
  PorousMediaBase(const InputParameters & parameters);
  virtual ~PorousMediaBase();
  

protected:
  void initStatefulProperties(unsigned n_points);
  virtual void computeProperties();

private:

  Real _input_default_diffusivity;
  Real _diffusivity_scaling_factor;

  Real _input_porosity;
  Real _input_bulk_density;
  Real _molecular_weight;
  Real _grain_density;

  Real _input_reactive_area;
  Real _rate_scaling_factor;
  

  Real _power_exponent;
  
  //some intrinsic material properties
  MaterialProperty<Real> & _diffusivity_of_O2;
  MaterialProperty<Real> & _diffusivity_of_N2;
  MaterialProperty<Real> & _diffusivity_of_CO;
  MaterialProperty<Real> & _diffusivity_of_CO2;
  MaterialProperty<Real> & _diffusivity_of_H2O;
  MaterialProperty<Real> & _diffusivity_of_H2;
  MaterialProperty<Real> & _diffusivity_of_NH3;
  MaterialProperty<Real> & _diffusivity_of_He;
  MaterialProperty<Real> & _diffusivity_of_Ar;
  

  MaterialProperty<Real> & _porosity;
  MaterialProperty<Real> & _porosity_old;
  
  MaterialProperty<Real> & _bulk_density;
  MaterialProperty<Real> & _bulk_density_old;

  MaterialProperty<Real> & _conversion_factor;
  

  MaterialProperty<Real> & _SA;         //reactive surface area
  MaterialProperty<Real> & _k_eff;      //reation rate constant
  MaterialProperty<Real> & _CO_to_CO2_ratio;  //reation product ratio
  
  

  const bool _has_temp; //coupled to temperature
  const VariableValue & _temperature;

  const bool _has_pressure; //coupled to pressure
  const VariableValue & _pressure;

  
  /// Gas constant, 8.314 J/mol/K
  Real _gas_const;
  /// System temperature
  Real _sys_temp;
  /// Actual system pressure
  Real _sys_pressure;

  
  /// Coupled gas mixture species concentrations
  std::vector<std::string> nl_vars;
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  std::vector<Real> _molecular_weights;
  std::vector<Real> _e_kb;
  std::vector<Real> _sigma;
  std::vector<Real> _dipole_moment;
  std::vector<Real> _polarizability;

  
   
 

#if LIBMESH_HAVE_PETSC
  Mat mat;
  Vec rhs;
  Vec sol;
  PetscScalar *values;
  PetscInt    *rows;
  PetscInt    *cols;
#endif

};

#endif //POROUSMEDIABASE_H
