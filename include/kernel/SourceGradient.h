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

#ifndef SOURCEGRADIENT_H
#define SOURCEGRADIENT_H

#include "Kernel.h"

class SourceGradient;

template<>
InputParameters validParams<SourceGradient>();


class SourceGradient : public Kernel
{
public:
  SourceGradient(const InputParameters & parameters);


  Real _coeff;

  const VariableValue & _coupled;

  const VariableValue & _coupled_old;

  const VariableGradient & _coupled_grad;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* SOURCEGRADIENT_H */
