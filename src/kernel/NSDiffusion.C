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

#include "NSDiffusion.h"


template<>
InputParameters validParams<NSDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<int>("direct", "dimensional component");
  params.addParam<MaterialPropertyName>("Reyold","Re",""); 

  return params;
}

NSDiffusion::NSDiffusion(const InputParameters & parameters) :
   Diffusion(parameters),
   _dim (getParam<int>("direct")),
   _Re(getMaterialProperty<Real>("Reyold"))

{}

Real
NSDiffusion::computeQpResidual()
{
	return _Re[_qp] * _grad_test[_i][_qp](_dim) * _grad_u[_qp](_dim);
}

Real
NSDiffusion::computeQpJacobian()
{
   return _Re[_qp] * _grad_test[_i][_qp](_dim) * _grad_phi[_j][_qp](_dim);
}

Real
NSDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}






