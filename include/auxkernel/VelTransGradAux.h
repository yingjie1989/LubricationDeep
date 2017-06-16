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

#ifndef VELTRANSGRADAUX_H
#define VELTRANSGRADAUX_H

#include "AuxKernel.h"

class VelTransGradAux;

template<>
InputParameters validParams<VelTransGradAux>();

class VelTransGradAux : public AuxKernel
{
public:

 VelTransGradAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableGradient & _coupled_grad;
  const VariableValue & _field;

  int _component;
  Real _value;

};


#endif //VELTRANSGRADAUX_H
