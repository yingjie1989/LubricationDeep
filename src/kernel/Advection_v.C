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

#include "Advection_v.h"


template<>
InputParameters validParams<Advection_v>()
{
  InputParameters params = validParams<Diffusion>();
//  params.addRequiredParam<std::string>("mat_name", "material property");
  params.addCoupledVar("coupled_h", "Coupled variable for thickness h");
  params.addCoupledVar("coupled_c", "Coupled variable for concentration c");
  params.addParam<int>("direct", "dimensional component");
  params.addParam<bool>("square", "indicator to couple h");
  params.addParam<MaterialPropertyName>("mu_param","mu","");
  params.addParam<MaterialPropertyName>("xi_param","xi","");
  params.addParam<MaterialPropertyName>("f_param","fmat","");
  params.addParam<Real>("sigma0_param","constant");
  params.addParam<Real>("power","power of surface tension"); 
  params.addParam<Real>("coeff","coefficient"); 

  return params;
}


Advection_v::Advection_v(const InputParameters & parameters) :
   Diffusion(parameters),
  // _mat_name(getParam<std::string>("mat_name")),
  // _mu(getMaterialPropertyByName<Real>(_mat_name + "mu_param")),
   _square(getParam<bool>("square")),
   _dim (getParam<int>("direct")),
   _coeff(getParam<Real>("coeff")),
   _power(getParam<Real>("power")),
   _coupled_h(coupledValue("coupled_h")),
   _coupled_c(coupledValue("coupled_c")),
   _coupled_gradc(coupledGradient("coupled_c")),
   _mu(getMaterialProperty<Real>("mu_param")),
   _xi(getMaterialProperty<Real>("xi_param")),
   _fmat(getMaterialProperty<Real>("f_param")),
  _sigma0(getParam<Real>("sigma0_param"))


{}

Real
Advection_v::computeQpResidual()
{

	if(_square)
             return -_coeff*_mu[_qp]*_test[_i][_qp]*(_coupled_h[_qp] - _fmat[_qp] )*(_coupled_h[_qp] - _fmat[_qp])*_coupled_gradc[_qp](_dim);
	else{
	    if(_coupled_c[_qp] < 1)
             	return - _power*pow(1-_coupled_c[_qp],_power-1)*_coeff*_xi[_qp]*_mu[_qp]*_test[_i][_qp]*(_coupled_h[_qp] - _fmat[_qp])*_coupled_gradc[_qp](_dim);
	     else 
	    	return 0.0;
	}
	return 0.0;
}

Real
Advection_v::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian


   return 0.0;
}

Real
Advection_v::computeQpOffDiagJacobian(unsigned int jvar)
{


 
  // get the coupled variable jvar is referring to

    int _coupled_h_aux = coupled("coupled_h");//check if this is the right component
    int _coupled_c_aux = coupled("coupled_c");

    if(_square){

        if (jvar == _coupled_h_aux)
                return -2*_coeff*_mu[_qp]*(_test[_i][_qp]*(_coupled_h[_qp] - _fmat[_qp] )*_phi[_j][_qp])*_coupled_gradc[_qp](_dim);
        if (jvar == _coupled_c_aux)
                return -_coeff*_mu[_qp]*_test[_i][_qp]*(_coupled_h[_qp]-_fmat[_qp])*(_coupled_h[_qp]-_fmat[_qp])*_grad_phi[_j][_qp](_dim);

	}

	else{
            if (jvar == _coupled_h_aux){
             if(_coupled_c[_qp] < 1)
                return -_power*pow(1-_coupled_c[_qp],_power-1)*_coeff*_xi[_qp]*_mu[_qp]*(_test[_i][_qp]*_phi[_j][_qp])*_coupled_gradc[_qp](_dim);
	    else 
		return 0.0;
	   }
          if (jvar == _coupled_c_aux){
	    if(_coupled_c[_qp] < 1)
                return -_coeff*_xi[_qp]*_mu[_qp]*_test[_i][_qp]*(_coupled_h[_qp]-_fmat[_qp])*( _grad_phi[_j][_qp](_dim)*_power*pow(1-_coupled_c[_qp],_power-1) - _coupled_gradc[_qp](_dim)*_power*(_power-1)*pow(1-_coupled_c[_qp],std::max(0.0,_power-2))*_phi[_j][_qp] );
	    else 
		return 0.0;
	  }
        }
	return 0.0;
}






