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

#ifndef ADVECTION_H_H
#define ADVECTION_H_H

#include "Diffusion.h"

class Advection_h;

template<>
InputParameters validParams<Advection_h>();


class Advection_h : public Diffusion
{
public:
  Advection_h(const InputParameters & parameters);


  //std::string _mat_name, _diffusion;

  const MaterialProperty<Real> & _fmat;
  const VariableValue & _coupled_x;
  const VariableValue & _coupled_y;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* ADVECTION_H_H */
