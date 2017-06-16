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

#include "NSSelf.h"


template<>
InputParameters validParams<NSSelf>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<int>("direct", "dimensional component");
  return params;
}


NSSelf::NSSelf(const InputParameters & parameters) :
   Diffusion(parameters),
  //   _Diffusivity(getMaterialProperty<Real>("Diffusivity"))
     _dim (getParam<int>("direct"))
{}

Real
NSSelf::computeQpResidual()
{
 // return _Diffusivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

    return _test[_i][_qp] * _u[_qp] * _grad_u[_qp](_dim);
}

Real
NSSelf::computeQpJacobian()
{
//  return _Diffusivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
 
    return _test[_i][_qp] * ( _u[_qp] * _grad_phi[_j][_qp](_dim) + _grad_u[_qp](_dim) * _phi[_j][_qp] ); 
}

Real
NSSelf::computeQpOffDiagJacobian(unsigned int jvar)
{
        return 0.0;
}

