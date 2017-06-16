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

#ifndef NSCOUPLED_H
#define NSCOUPLED_H

#include "Diffusion.h"

class NSCoupled;

template<>
InputParameters validParams<NSCoupled>();


class NSCoupled : public Diffusion
{
public:
  NSCoupled(const InputParameters & parameters);


  //std::string _mat_name;
  
  int  _dim;
  const VariableValue & _coupled_v;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* NSCOUPLED_H */
