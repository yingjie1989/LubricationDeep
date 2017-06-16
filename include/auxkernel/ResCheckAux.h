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

#ifndef RESCHECKAUX_H
#define RESCHECKAUX_H

#include "AuxKernel.h"

class ResCheckAux;

template<>
InputParameters validParams<ResCheckAux>();

class ResCheckAux : public AuxKernel
{
public:

 ResCheckAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled_u;
  const VariableValue & _coupled_v;
  const VariableValue & _coupled_w;

  const VariableGradient & _coupled_gradu;
  const VariableGradient & _coupled_gradv;
  const VariableGradient & _coupled_gradw;

  Real _value;

};


#endif //RESCHECKAUX_H
