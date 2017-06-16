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

#include "AdvectionSUPG.h"

template<>
InputParameters validParams<AdvectionSUPG>()
{
  InputParameters params = validParams<Diffusion>();;
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addCoupledVar("coupled_x", 0, "Coupled variable for component x");
  params.addCoupledVar("coupled_y", 0, "Coupled variable for component y");
  params.addCoupledVar("coupled_z", 0, "Coupled variable for component z");
  params.addParam<Real>("coeff_param",1,"additional viscous term");
  params.addParam<int>("dim",1,"dimension");
  params.addParam<Real>("tol_param",0,"tolerance");
  params.addParam<MaterialPropertyName>("mesh_size","hele","meshSize");  
  params.addParam<MaterialPropertyName>("D_param","DParameter","");
  return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters) :
    Diffusion(parameters),
    _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second()),
    _coupled_x(coupledValue("coupled_x")),	
    _coupled_y(coupledValue("coupled_y")),
    _coupled_z(coupledValue("coupled_z")),
    _coupled_grad_x(coupledGradient("coupled_x")),
    _coupled_grad_y(coupledGradient("coupled_y")),
    _coeff(getParam<Real>("coeff_param")),
    _dim(getParam<int>("dim")),
    _tol(getParam<Real>("tol_param")),
    _hele(getMaterialProperty<Real>("mesh_size")), 
    _Diff(getMaterialProperty<Real>("D_param"))

{}

Real
AdvectionSUPG::computeQpResidual()
{
 
   Real _uN = std::sqrt(_coupled_x[_qp]*_coupled_x[_qp] + _coupled_y[_qp]*_coupled_y[_qp]);
  
   Real _tau = 0.0;
   Real _psi = 0.0;
   Real _coeff1 = 1.0;
   Real _coeff2 = 0.0;

   if(_Diff[_qp]>0){
      //  Real _alpha = _uN*_hele[_qp]/_Diff[_qp]*0.5;
      //  if(_alpha < 1e-3) _psi = _alpha/3.0;
      //  else _psi = 1.0/std::tanh(_alpha) - 1.0/_alpha;
      //  _tau = _hele[_qp]*_psi*0.5;

	 _tau = 0.5*_hele[_qp]*_uN/sqrt( 4.0*_Diff[_qp]*_Diff[_qp]/(_hele[_qp]*_hele[_qp]) + _uN*_uN );

   }else{
        _tau = _hele[_qp]*0.5;
	_coeff1 = 0.0;	
	_coeff2 = 1.0;

    }


    Real signu = 0.0;
    
    if(_coupled_x[_qp] > _tol) signu =  1.0;
    if(_coupled_x[_qp] < -_tol) signu = -1.0;
  
    if(_dim == 2)
	return _coeff * _tau * ( _grad_test[_i][_qp](0) * (_coupled_x[_qp]/_uN) + _grad_test[_i][_qp](1) * ( _coupled_y[_qp]/_uN)  ) * ( _grad_u[_qp](0)*_coupled_x[_qp] + _grad_u[_qp](1)*_coupled_y[_qp]  - _coeff2 * _coupled_z[_qp] + _coeff1 * _u[_qp]*(_coupled_grad_x[_qp](0) + _coupled_grad_y[_qp](1) ) - _Diff[_qp] * ( _second_u[_qp](0,0) + _second_u[_qp](1,1) )  ); 
 
    return _coeff * _tau*( _grad_test[_i][_qp](0)*signu )*( _grad_u[_qp](0)*_coupled_x[_qp] + _grad_u[_qp](1)*_coupled_y[_qp]  - _coeff2 * _coupled_z[_qp] + _coeff1 * _u[_qp]*(_coupled_grad_x[_qp](0) + _coupled_grad_y[_qp](1) ) - _Diff[_qp] * ( _second_u[_qp](0,0) + _second_u[_qp](1,1) ) ); 


}

Real
AdvectionSUPG::computeQpJacobian()
{

   Real _uN = std::sqrt(_coupled_x[_qp]*_coupled_x[_qp] + _coupled_y[_qp]*_coupled_y[_qp]);

   Real _tau = 0.0;
   Real _psi = 0.0;
   Real _coeff1 = 1.0;	  
   Real _coeff2 = 0.0;


   if(_Diff[_qp]>0){
      //  Real _alpha = _uN*_hele[_qp]/_Diff[_qp]*0.5;
      //  if(_alpha < 1e-3) _psi = _alpha/3.0;
      //  else _psi = 1.0/std::tanh(_alpha) - 1.0/_alpha;
      //  _tau = _hele[_qp]*_psi*0.5;

	_tau = 0.5*_hele[_qp]*_uN/sqrt( 4.0*_Diff[_qp]*_Diff[_qp]/(_hele[_qp]*_hele[_qp]) + _uN*_uN );

	//_tau = _hele[_qp]*0.5;
   }else{
        _tau = _hele[_qp]*0.5;
	_coeff1 = 0.0;
	_coeff2 = 1.0;

   }

	 Real signu = 0.0;

	if(_coupled_x[_qp] > _tol) signu =  1.0;
    	if(_coupled_x[_qp] < -_tol) signu = -1.0;

	if(_dim == 2)
		return _coeff * _tau * ( _grad_test[_i][_qp](0) * (_coupled_x[_qp]/_uN) + _grad_test[_i][_qp](1) * ( _coupled_y[_qp]/_uN)  ) * ( _grad_phi[_j][_qp](0)*(_coupled_x[_qp]) + _grad_phi[_j][_qp](1)*(_coupled_y[_qp]) + _coeff1 * _phi[_j][_qp]*(_coupled_grad_x[_qp](0) + _coupled_grad_y[_qp](1) ) - _Diff[_qp] * ( _second_phi[_j][_qp](0,0) + _second_phi[_j][_qp](1,1) ) );
	
	
	return _coeff * _tau*( _grad_test[_i][_qp](0)*signu )*( _grad_phi[_j][_qp](0)*(_coupled_x[_qp]) + _grad_phi[_j][_qp](1)*(_coupled_y[_qp]) + _coeff1 * _phi[_j][_qp]*(_coupled_grad_x[_qp](0) + _coupled_grad_y[_qp](1) ) - _Diff[_qp] * ( _second_phi[_j][_qp](0,0) + _second_phi[_j][_qp](1,1) ) );

	

}


Real
AdvectionSUPG::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}

