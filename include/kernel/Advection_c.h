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

#ifndef ADVECTION_C_H
#define ADVECTION_C_H

#include "Diffusion.h"

class Advection_c;

template<>
InputParameters validParams<Advection_c>();


class Advection_c : public Diffusion
{
public:
  Advection_c(const InputParameters & parameters);


  //std::string _mat_name, _diffusion;

  int _dim;
  const VariableValue & _coupled_x;
  const VariableGradient & _coupled_grad_x;
  const VariableValue & _coupled_y;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* ADVECTION_C_H */
