/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef InNSSUPGMOMENTUM_H
#define InNSSUPGMOMENTUM_H

#include "Kernel.h"

// Forward Declarations
class InNSSUPGMomentum;

template<>
InputParameters validParams<InNSSUPGMomentum>();

/**
 * This class computes momentum equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation.
 */
class InNSSUPGMomentum : public Kernel
{
public:
  InNSSUPGMomentum(const InputParameters & parameters);

  virtual ~InNSSUPGMomentum(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

   const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _p;

  const VariableValue & _u_vel_mesh;
  const VariableValue & _v_vel_mesh;
  const VariableValue & _w_vel_mesh;


  const VariableValue & _u_vel_old;
  const VariableValue & _v_vel_old;
  const VariableValue & _w_vel_old;


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
  Real _DeltaT;
  Real _mu;
  Real _visco_coeff;
  Real _grav;
  Real _rho;
  std::vector<Real> _meshsize;
  bool _ifOld;
  // Parameters
  unsigned _component;
  Real _coeff;
  Real _coeff_param;
  unsigned _vert;
};


#endif // INSSUPGMOMENTUM_H
