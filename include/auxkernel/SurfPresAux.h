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

#ifndef SURFPRESAUX_H
#define SURFPRESAUX_H

#include "AuxKernel.h"
#include "Assembly.h"


class SurfPresAux;

template<>
InputParameters validParams<SurfPresAux>();

class SurfPresAux : public AuxKernel
{
public:

 SurfPresAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled_c; 
  const VariableGradient & _grad_h;
  const VariableSecond & _h_second;
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  unsigned _dim;
  Real _gamma0;
  Real _mu;

};


#endif //SURFPRESAUX_H
