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

#ifndef NSDIFFUSION_H
#define NSDIFFUSION_H

#include "Diffusion.h"

class NSDiffusion;

template<>
InputParameters validParams<NSDiffusion>();


class NSDiffusion : public Diffusion
{
public:
  NSDiffusion(const InputParameters & parameters);


  //std::string _mat_name;
  
  int  _dim;
  const MaterialProperty<Real> & _Re;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* NSDIFFUSION_H */
