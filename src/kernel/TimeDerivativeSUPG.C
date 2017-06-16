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

#include "TimeDerivativeSUPG.h"
#include "Assembly.h"

#include "libmesh/quadrature.h"

template<>
InputParameters validParams<TimeDerivativeSUPG>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addCoupledVar("coupled_x", 0, "Coupled variable for component x");
  params.addCoupledVar("coupled_y", 0, "Coupled variable for component y");
  params.addParam<Real>("tol_param",0,"tolerance");
  params.addParam<int>("dim",1,"dimension");
  params.addParam<Real>("coeff_param",1,"addtional viscous term");
  params.addParam<MaterialPropertyName>("mesh_size","hele","meshSize");
  params.addParam<MaterialPropertyName>("D_param",0,"");
  return params;
}

TimeDerivativeSUPG::TimeDerivativeSUPG(const InputParameters & parameters) :
    TimeKernel(parameters),
    _coupled_x(coupledValue("coupled_x")),	
    _coupled_y(coupledValue("coupled_y")),
    _coupled_old_x(coupledValueOld("coupled_x")),
    _coupled_old_y(coupledValueOld("coupled_y")),
    _coupled_old_grad_x(coupledGradientOld("coupled_x")),
    _coupled_old_grad_y(coupledGradientOld("coupled_y")), 
    _coeff(getParam<Real>("coeff_param")),
    _dim(getParam<int>("dim")),
    _tol(getParam<Real>("tol_param")),
    _hele(getMaterialProperty<Real>("mesh_size")),
    _Diff(getMaterialProperty<Real>("D_param")),
    _lumping(getParam<bool>("lumping"))
{
}

Real
TimeDerivativeSUPG::computeQpResidual()
{
    
    Real _uN = std::sqrt(_coupled_x[_qp]*_coupled_x[_qp] + _coupled_y[_qp]*_coupled_y[_qp]);
  
    Real _tau = 0.0;
    Real _psi = 0.0; 

    if(_Diff[_qp]>0){       
    	//Real _alpha = _uN*_hele[_qp]/_Diff[_qp]*0.5;
        //if(_alpha < 1e-3) _psi = _alpha/3.0;
        //else _psi = 1.0/std::tanh(_alpha) - 1.0/_alpha;
	
	//_tau = _hele[_qp]*_psi*0.5;

	 _tau = 0.5*_hele[_qp]*_uN/sqrt( 4.0*_Diff[_qp]*_Diff[_qp]/(_hele[_qp]*_hele[_qp]) + _uN*_uN );


   }else{
	_tau = _hele[_qp]*0.5;
    }
		
    Real signu = 0.0;

    if(_coupled_x[_qp] > _tol) signu =  1.0;
    if(_coupled_x[_qp] < -_tol) signu = -1.0;

    if(_dim == 2)
	return _coeff * _tau * ( _grad_test[_i][_qp](0) * (_coupled_x[_qp]/_uN) + _grad_test[_i][_qp](1) * ( _coupled_y[_qp]/_uN)  ) * _u_dot[_qp]; 

    return _coeff * _tau*( _grad_test[_i][_qp](0)*signu )*_u_dot[_qp];
} 

Real
TimeDerivativeSUPG::computeQpJacobian()
{

    Real _uN = std::sqrt(_coupled_x[_qp]*_coupled_x[_qp] + _coupled_y[_qp]*_coupled_y[_qp]);
         
    Real _tau = 0.0;
    Real _psi = 0.0;   
	
    if(_Diff[_qp]>0){
       // Real _alpha = _uN*_hele[_qp]/_Diff[_qp]*0.5;
       // if(_alpha < 1e-3) _psi = _alpha/3.0;
       // else _psi = 1.0/std::tanh(_alpha) - 1.0/_alpha;
       // _tau = _hele[_qp]*_psi*0.5;

      _tau = 0.5*_hele[_qp]*_uN/sqrt( 4.0*_Diff[_qp]*_Diff[_qp]/(_hele[_qp]*_hele[_qp]) + _uN*_uN );

   }else{
        _tau = _hele[_qp]*0.5;
  }

    Real signu = 0.0;

    if(_coupled_x[_qp] > _tol) signu =  1.0;
    if(_coupled_x[_qp] < -_tol) signu = -1.0;

    if(_dim == 2)
        return _coeff * _tau * ( _grad_test[_i][_qp](0) * (_coupled_x[_qp]/_uN) + _grad_test[_i][_qp](1) * ( _coupled_y[_qp]/_uN)  ) * _phi[_j][_qp]*_du_dot_du[_qp]; 


    return _coeff * _tau*( _grad_test[_i][_qp](0)*signu )*_phi[_j][_qp]*_du_dot_du[_qp];


}


Real
TimeDerivativeSUPG::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}

