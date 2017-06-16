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

#ifndef DEFORMAUX_H
#define DEFORMAUX_H

#include "AuxKernel.h"

class DeformAux;

template<>
InputParameters validParams<DeformAux>();

class DeformAux : public AuxKernel
{
public:

 DeformAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled_val;

  Real _value;
};


#endif //DEFORMAUX_H
