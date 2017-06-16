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

#ifndef DIFFUSIONCOUPLED_H
#define DIFFUSIONCOUPLED_H

#include "Diffusion.h"

class DiffusionCoupled;

template<>
InputParameters validParams<DiffusionCoupled>();


class DiffusionCoupled : public Diffusion
{
public:
  DiffusionCoupled(const InputParameters & parameters);


  //std::string _mat_name, _diffusion;

  const VariableValue & _coupled_c;
  const VariableGradient & _coupled_grad_c;
  const VariableGradient & _coupled_grad_h;
  const MaterialProperty<Real> & _xi;
  unsigned _vert;
  Real _power;
  Real _sigma0;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* DIFFUSIONCOUPLED_H */
