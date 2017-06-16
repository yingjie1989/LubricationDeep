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

#include "DiffusionCoupled.h"


template<>
InputParameters validParams<DiffusionCoupled>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<MaterialPropertyName>("xi_param", "xi", "surface tension parameter");
  params.addCoupledVar("coupled_c", 0.0, "Coupled variable, concentration");
  params.addCoupledVar("coupled_h",0.0, "Coupled variable, height");
  params.addParam<unsigned>("dim",2,"dimension");
  params.addParam<Real>("power","suracepower");
  params.addParam<Real>("sigma0_param","constant");
  return params;
}


DiffusionCoupled::DiffusionCoupled(const InputParameters & parameters) :
   Diffusion(parameters),

     _coupled_c(coupledValue("coupled_c")),
     _coupled_grad_c(coupledGradient("coupled_c")),
     _coupled_grad_h(coupledGradient("coupled_h")),
     _xi(getMaterialProperty<Real>("xi_param")),
     _vert(getParam<unsigned>("dim")),
     _power(getParam<Real>("power")),
     _sigma0(getParam<Real>("sigma0_param"))

{}

Real
DiffusionCoupled::computeQpResidual()
{

  Real value =  _coupled_grad_h[_qp](0) * (_grad_test[_i][_qp](0)*(_sigma0 + _xi[_qp]*pow((1.0-_coupled_c[_qp]),_power))-_test[_i][_qp]*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*_coupled_grad_c[_qp](0));

  if(_vert > 1)
         value += _coupled_grad_h[_qp](1) * (_grad_test[_i][_qp](1)*(_sigma0 + _xi[_qp]*pow((1.0-_coupled_c[_qp]),_power))-_test[_i][_qp]*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*_coupled_grad_c[_qp](1));

	return value;
}

Real
DiffusionCoupled::computeQpJacobian()
{
  return 0.0;
}


Real
DiffusionCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
    int _coupled_c_aux = coupled("coupled_c");
    int _coupled_h_aux = coupled("coupled_h");
    if (jvar==_coupled_c_aux){

	Real value =  _coupled_grad_h[_qp](0) * ( _grad_test[_i][_qp](0)*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*(-_phi[_j][_qp])  
	- _test[_i][_qp]*_xi[_qp]*_power*( (_power-1)*pow((1.0-_coupled_c[_qp]),std::max(_power-2,0.0))*_coupled_grad_c[_qp](0)*(-_phi[_j][_qp]) + pow((1.0-_coupled_c[_qp]),_power-1)* _grad_phi[_j][_qp](0)  ) );
	if(_vert > 1)

	value += _coupled_grad_h[_qp](1) * ( _grad_test[_i][_qp](1)*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*(-_phi[_j][_qp])  - _test[_i][_qp]*_xi[_qp]*_power*( (_power-1)*pow((1.0-_coupled_c[_qp]),std::max(_power-2,0.0))*_coupled_grad_c[_qp](1)*(-_phi[_j][_qp]) + pow((1.0-_coupled_c[_qp]),_power-1)* _grad_phi[_j][_qp](1)  ) );

	return value;
}
    if (jvar==_coupled_h_aux)
{
        Real value = _grad_phi[_j][_qp](0) * (_grad_test[_i][_qp](0)*(_sigma0 + _xi[_qp]*pow((1.0-_coupled_c[_qp]),_power))-_test[_i][_qp]*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*_coupled_grad_c[_qp](0));
	
	if(_vert > 1)

	     value += _grad_phi[_j][_qp](1) * (_grad_test[_i][_qp](1)*(_sigma0 + _xi[_qp]*pow((1.0-_coupled_c[_qp]),_power))-_test[_i][_qp]*_xi[_qp]*_power*pow((1.0-_coupled_c[_qp]),_power-1)*_coupled_grad_c[_qp](1));

	return value;
    }


  return 0.0;
}










