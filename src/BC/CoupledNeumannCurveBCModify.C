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

#include "CoupledNeumannCurveBCModify.h"

template<>
InputParameters validParams<CoupledNeumannCurveBCModify>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol", "indicate tolerance");
        params.addParam<bool>("old", false,  "indicate dimension");
        params.addParam<unsigned>("direct", "indicate dimension");
        params.addParam<unsigned>("component", "indicate component");
        params.addCoupledVar("coupled_p",0.0, "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_c", "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_h", "Flux Value at the Boundary");

        return params;
}

CoupledNeumannCurveBCModify::CoupledNeumannCurveBCModify(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _ifOld(getParam<bool>("old")),
    _tol(getParam<Real>("tol")),   
    _dim (getParam<unsigned>("direct")),
    _component(getParam<unsigned>("component")), 
    _some_var_val(coupledValue("coupled_c")),
    _p_coupled(coupledValue("coupled_p")),
    _p_coupled_old(coupledValueOld("coupled_p")),
    _c_grad(coupledGradient("coupled_c")),
    _grad_h(coupledGradient("coupled_h")), 
    _h_second(coupledSecond("coupled_h"))
    

{}

Real
CoupledNeumannCurveBCModify::computeQpResidual()
{

   Real norm2d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) );
   RealVectorValue normal2d(-_grad_h[_qp](0)/norm2d, 1/norm2d, 0);
   RealVectorValue tau2d(1/norm2d,_grad_h[_qp](0)/norm2d,0);

   Real norm3d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue normal3d(-_grad_h[_qp](0)/norm3d, -_grad_h[_qp](1)/norm3d, 1/norm3d);
   Real normp = sqrt( _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue tau3d(_grad_h[_qp](0)/norm3d/normp, _grad_h[_qp](1)/norm3d/normp, normp/norm3d);

   Real sigma_val = 0.0;
   Real _sigma_n = 0.0;
   Real _sigma_tx = 0.0;
   Real _sigma_ty = 0.0;

   if(_some_var_val[_qp] < _tol )
	return 0.0;


   if(_dim == 1){
	_sigma_tx = -_Re*( _c_grad[_qp](0) * tau2d(0) + _c_grad[_qp](1) * tau2d(1) );

	Real gamma = _Re*( 1 - _some_var_val[_qp] );
        Real Rinv = _h_second[_qp](0,0);//(norm2d*norm2d*norm2d);    

        _sigma_n  = gamma * Rinv; 	

	if(_component == 0 )
		sigma_val = _sigma_n * normal2d(0) + _sigma_tx * tau2d(0);
	else if(_component == 1)
		sigma_val = _sigma_n * normal2d(1) + _sigma_tx * tau2d(1);
	else
		sigma_val = 0.0;

   }else if(_dim == 2){ 
        
        Real _sigma_ts = -_Re*( _c_grad[_qp](0) * tau3d(0) + _c_grad[_qp](1) * tau3d(1) + _c_grad[_qp](2) * tau3d(2) );
	
	_sigma_tx = _sigma_ts * _grad_h[_qp](0)/normp;
	_sigma_ty = _sigma_ts * _grad_h[_qp](1)/normp;

        Real gamma = _Re * (1 - _some_var_val[_qp]);
        Real Rinv = 0.5*( (1 + _grad_h[_qp](0)*_grad_h[_qp](0))*_h_second[_qp](1,1) + (1 + _grad_h[_qp](1)*_grad_h[_qp](1))*_h_second[_qp](0,0) - 2*_grad_h[_qp](0)*_grad_h[_qp](1)*_h_second[_qp](0,1) )/(norm3d*norm3d*norm3d);

        _sigma_n  = gamma * Rinv;    
  
	if(_component == 0)
		sigma_val = _sigma_tx; 
	else if(_component == 1)
		sigma_val = _sigma_ty;
	else if(_component == 2)
		sigma_val = _sigma_n;
	else
		sigma_val =  0.0;
   }
	return - _test[_i][_qp] * sigma_val;
}
