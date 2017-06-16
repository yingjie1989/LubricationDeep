/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef InNSMOMENTUM_H
#define InNSMOMENTUM_H

#include "Kernel.h"

// Forward Declarations
class InNSMomentum;

template<>
InputParameters validParams<InNSMomentum>();

/**
 * This class computes momentum equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation.
 */
class InNSMomentum : public Kernel
{
public:
  InNSMomentum(const InputParameters & parameters);

  virtual ~InNSMomentum(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _p;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;
  const VariableGradient & _grad_p;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _p_var_number;

  // Material properties
  // MaterialProperty<Real> & _dynamic_viscosity;
  Real _mu;
  Real _rho;
  Real _gravity;
  Real _coeff;

  // Parameters
  unsigned _component;
  unsigned _vert;
  bool _pressure;
};


#endif // INSMOMENTUM_H
