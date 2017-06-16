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

#ifndef ADVECTION_V_H
#define ADVECTION_V_H

#include "Diffusion.h"

class Advection_v;

template<>
InputParameters validParams<Advection_v>();


class Advection_v : public Diffusion
{
public:
  Advection_v(const InputParameters & parameters);


  //std::string _mat_name;
  
  bool _square;
  int  _dim;
  Real  _coeff;
  Real _power;
  const VariableValue & _coupled_h;
  const VariableValue & _coupled_c;
  const VariableGradient & _coupled_gradc;
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _xi;
  const MaterialProperty<Real> & _fmat;
  Real _sigma0;


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* ADVECTION_v_H */
