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

#ifndef INCOMPMODIFY_H
#define INCOMPMODIFY_H

#include "Kernel.h"

class InCompModify;

template<>
InputParameters validParams<InCompModify>();


class InCompModify : public Kernel
{
public:
  InCompModify(const InputParameters & parameters);
  int  _dim;
  Real _mu;
  const VariableGradient & _coupled_gradv;
  const VariableGradient & _coupled_gradu;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* INCOMPMODIFY */
