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

#ifndef NSSELF_H
#define NSSELF_H

#include "Diffusion.h"

class NSSelf;

template<>
InputParameters validParams<NSSelf>();


class NSSelf : public Diffusion
{
public:
  NSSelf(const InputParameters & parameters);

//  const MaterialProperty<Real> & _Diffusivity;
  int _dim;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);

};


#endif /* NSSELF_H */
