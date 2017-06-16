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

#ifndef INCOMPRESSVELOCITY_H
#define INCOMPRESSVELOCITY_H

#include "Diffusion.h"

class InCompressVelocity;

template<>
InputParameters validParams<InCompressVelocity>();


class InCompressVelocity : public Diffusion
{
public:
  InCompressVelocity(const InputParameters & parameters);

  
  int  _dim;
  //const VariableValue & _coupled_v;
  const VariableGradient & _coupled_gradv;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* INCOMPRESSVELOCITY */
