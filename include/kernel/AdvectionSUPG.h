/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ADVECTIONSUPG_H
#define ADVECTIONSUPG_H

#include "Diffusion.h"

// Forward Declaration
class AdvectionSUPG;

template<>
InputParameters validParams<AdvectionSUPG>();

class AdvectionSUPG : public Diffusion
{
public:
  AdvectionSUPG(const InputParameters & parameters);

  //virtual void computeJacobian();
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

  const VariableValue & _coupled_x;
  const VariableValue & _coupled_y;
  const VariableValue & _coupled_z;
  const VariableGradient & _coupled_grad_x;
  const VariableGradient & _coupled_grad_y;

  Real _coeff;
  int _dim;
  Real _tol;
  const MaterialProperty<Real> & _hele;
  const MaterialProperty<Real> & _Diff;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
};

#endif //ADVECTIONSUPG_H
