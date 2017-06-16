/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef InNSSUPGTIMEDERIVATIVE_H
#define InNSSUPGTIMEDERIVATIVE_H

#include "TimeKernel.h"

// Forward Declarations
class InNSSUPGTimeDerivative;

template<>
InputParameters validParams<InNSSUPGTimeDerivative>();

/**
 * This class computes momentum equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation.
 */
class InNSSUPGTimeDerivative : public TimeKernel
{
public:
  InNSSUPGTimeDerivative(const InputParameters & parameters);

  virtual ~InNSSUPGTimeDerivative(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
 
   const VariableValue & _u_vel_old;
  const VariableValue & _v_vel_old;
  const VariableValue & _w_vel_old;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;

  const VariableValue & _u_vel_mesh;
  const VariableValue & _v_vel_mesh;
  const VariableValue & _w_vel_mesh;


  // Material properties
  // MaterialProperty<Real> & _dynamic_viscosity;
  Real _DeltaT;
  Real _mu;
  Real _visco_coeff;
  Real _rho;
  std::vector<Real> _meshsize;
  bool _ifOld;
  // Parameters
  unsigned _component;
  unsigned _vert;
};


#endif // INNSSUPGTIMEDERIVATIVE_H
