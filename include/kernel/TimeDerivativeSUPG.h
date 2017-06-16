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

#ifndef TIMEDERIVATIVESUPG_H
#define TIMEDERIVATIVESUPG_H

#include "TimeKernel.h"

// Forward Declaration
class TimeDerivativeSUPG;

template<>
InputParameters validParams<TimeDerivativeSUPG>();

class TimeDerivativeSUPG : public TimeKernel
{
public:
  TimeDerivativeSUPG(const InputParameters & parameters);

  //virtual void computeJacobian();
  const VariableValue & _coupled_x;
  const VariableValue & _coupled_y;
  const VariableValue & _coupled_old_x;
  const VariableValue & _coupled_old_y;
  const VariableGradient & _coupled_old_grad_x;
  const VariableGradient & _coupled_old_grad_y;

  Real _coeff;
  int _dim;
  Real _tol;
  const MaterialProperty<Real> & _hele;
  const MaterialProperty<Real> & _Diff;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  bool _lumping;
};

#endif //TIMEDERIVATIVESUPG_H
