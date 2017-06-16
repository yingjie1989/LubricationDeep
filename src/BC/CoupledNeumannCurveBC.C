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

#include "CoupledNeumannCurveBC.h"

template<>
InputParameters validParams<CoupledNeumannCurveBC>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol", "indicate tolerance");
        params.addParam<bool>("old", false,  "indicate dimension");
        params.addParam<int>("direct", "indicate dimension");
        params.addCoupledVar("coupled_p",0.0, "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_c", "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_h", "Flux Value at the Boundary");

        return params;
}

CoupledNeumannCurveBC::CoupledNeumannCurveBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _ifOld(getParam<bool>("old")),
    _tol(getParam<Real>("tol")),   
    _dim (getParam<int>("direct")), 
    _some_var_val(coupledValue("coupled_c")),
    _p_coupled(coupledValue("coupled_p")),
    _p_coupled_old(coupledValueOld("coupled_p")),
    _c_grad(coupledGradient("coupled_c")),
    _h_grad(coupledGradient("coupled_h")), 
    _h_second(coupledSecond("coupled_h"))
    

{}

Real
CoupledNeumannCurveBC::computeQpResidual()
{

	Real _p_value = 0.0;
  if(_ifOld)
	_p_value = _p_coupled_old[_qp];
  else
	_p_value = _p_coupled[_qp];

  if(_some_var_val[_qp] > _tol){	
	Real tau = sqrt(1 + _h_grad[_qp](_dim) * _h_grad[_qp](_dim));
	Real gamma = 1 - _some_var_val[_qp];
	Real Rinv = _h_second[_qp](_dim,_dim)/(tau*tau*tau);
	
  	return - 0.5 * _test[_i][_qp]*_Re*( gamma*Rinv + _p_value );
   }
  return 0.0;
}
