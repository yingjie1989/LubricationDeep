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

#include "CoupledNeumannGradientBCModify.h"

template<>
InputParameters validParams<CoupledNeumannGradientBCModify>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol", "indicate tolerance");
        params.addParam<int>("dimension", "indicate dimension");
        params.addParam<int>("direct", "indicate direction");
        params.addCoupledVar("some_var",0, "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_c", "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_h", "Flux Value at the Boundary");

        return params;
}

CoupledNeumannGradientBCModify::CoupledNeumannGradientBCModify(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _tol(getParam<Real>("tol")),   
    _dim(getParam<int>("dimension")),
    _dirt(getParam<int>("direct")), 
    _some_var_val(coupledValue("some_var")),
    _c_grad(coupledGradient("coupled_c")),
    _h_grad(coupledGradient("coupled_h"))
{}

Real
CoupledNeumannGradientBCModify::computeQpResidual()
{
  if(_dim == 2){	
	Real tau = sqrt(1 + _h_grad[_qp](0) * _h_grad[_qp](0));
  	return -_test[_i][_qp]*_Re*( _c_grad[_qp](0) + _c_grad[_qp](1)*_h_grad[_qp](0) )/tau;
   
  }else{
	 Real norm2d = sqrt( 1 + _h_grad[_qp](0) * _h_grad[_qp](0) );
	 Real norm3d = sqrt( 1 + _h_grad[_qp](0) * _h_grad[_qp](0) + _h_grad[_qp](1) * _h_grad[_qp](1) );
   	 RealVectorValue normal3d(-_h_grad[_qp](0)/norm3d, -_h_grad[_qp](1)/norm3d, 1/norm3d);
   	 Real norm2dy = sqrt( 1 + _h_grad[_qp](1) * _h_grad[_qp](1) );
   	 RealVectorValue tau3dx(1/norm2d, 0, _h_grad[_qp](0)/norm2d);
   	 RealVectorValue tau3dy(0, 1/norm2dy, _h_grad[_qp](1)/norm2dy);
	
	if(_dirt == 0)
	  	return  - _test[_i][_qp] * _Re * ( _c_grad[_qp] * tau3dx );
        else if(_dirt == 1)
		return  - _test[_i][_qp] * _Re * ( _c_grad[_qp] * tau3dy );
	else 
		return 0.0;
  }
  return 0.0;
}
