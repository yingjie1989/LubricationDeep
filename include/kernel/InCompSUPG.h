/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef InCOMPSUPG_H
#define InCOMPSUPG_H

#include "Kernel.h"

// Forward Declarations
class InCompSUPG;

template<>
InputParameters validParams<InCompSUPG>();

/**
 * This class computes momentum equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation.
 */
class InCompSUPG : public Kernel
{
public:
  InCompSUPG(const InputParameters & parameters);

  virtual ~InCompSUPG(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  const MaterialProperty<Real> & _grav;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  const VariableValue & _u_vel_mesh;
  const VariableValue & _v_vel_mesh;
  const VariableValue & _w_vel_mesh;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  //Second Derivative

  const VariableSecond & _second_u_vel;
  const VariableSecond & _second_v_vel;
  const VariableSecond & _second_w_vel;
  const VariablePhiSecond & _second_phi;

 

  // TimeDerivative

  const VariableValue & _u_vel_dot;
  const VariableValue & _v_vel_dot;
  const VariableValue & _w_vel_dot;

  // TimeDerivative with respect to the coefficients

  const VariableValue & _u_vel_dot_du;
  const VariableValue & _v_vel_dot_du;
  const VariableValue & _w_vel_dot_du;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;

  // Material properties
  // MaterialProperty<Real> & _dynamic_viscosity;
  Real _DeltaT;
  std::vector<Real> _mu;
  std::vector<Real> _rho;
  std::vector<Real> _meshsize;

  // Parameters
  unsigned _component;
  Real _coeff_pressure;
  Real _coeff_param; 
  unsigned _vert;
};


#endif // InCOMPSUPG_H
