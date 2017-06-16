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

#include "CoupledNeumannGradientBC.h"

template<>
InputParameters validParams<CoupledNeumannGradientBC>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol", "indicate tolerance");
        params.addParam<int>("direct", "indicate dimension");
        params.addRequiredCoupledVar("some_var", "Flux Value at the Boundary");
        return params;
}

CoupledNeumannGradientBC::CoupledNeumannGradientBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _tol(getParam<Real>("tol")),   
    _dim (getParam<int>("direct")), 
    _some_var_val(coupledValue("some_var")),
    _some_var_grad(coupledGradient("some_var"))
{}

Real
CoupledNeumannGradientBC::computeQpResidual()
{
  if(_some_var_val[_qp] > _tol)	
  	return -_test[_i][_qp]*_Re*_some_var_grad[_qp](_dim);
  return 0.0;
}
