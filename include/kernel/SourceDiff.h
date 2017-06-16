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

#ifndef SOURCEDIFF_H
#define SOURCEDIFF_H

#include "Diffusion.h"

class SourceDiff;

template<>
InputParameters validParams<SourceDiff>();


class SourceDiff : public Diffusion
{
public:
  SourceDiff(const InputParameters & parameters);


  Real _component;
  Real _coeff;

  const VariableValue & _coupled;

  const VariableValue & _coupled_old;

  const VariableGradient & _coupled_grad;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* SOURCEDIFF_H */
