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

#ifndef INCOMPRESS_H
#define INCOMPRESS_H

#include "Kernel.h"

class InCompress;

template<>
InputParameters validParams<InCompress>();


class InCompress : public Kernel
{
public:
  InCompress(const InputParameters & parameters);
  int  _dim;
  const VariableGradient & _coupled_gradv;
  const VariableGradient & _coupled_gradu;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* INCOMPRESS */
