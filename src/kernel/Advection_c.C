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

#include "Advection_c.h"


template<>
InputParameters validParams<Advection_c>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<int>("dim", "dimension of problem");   
  params.addCoupledVar("coupled_x",0, "Coupled variable for component x");
  params.addCoupledVar("coupled_y",0, "Coupled variable for component y");
  //params.addRequiredParam<std::string>("mat_name", "Praft material name");

  return params;
}


Advection_c::Advection_c(const InputParameters & parameters) :
   Diffusion(parameters),
    //_mat_name(getParam<std::string>("mat_name")),
    _dim (getParam<int>("dim")),
    _coupled_x(coupledValue("coupled_x")),
    _coupled_grad_x(coupledGradient("coupled_x")),
    _coupled_y(coupledValue("coupled_y"))

{}

Real
Advection_c::computeQpResidual()
{
       if(_dim == 2)
           return -(_grad_test[_i][_qp](0)*_coupled_x[_qp] + _grad_test[_i][_qp](1)*_coupled_y[_qp])*_u[_qp];
       else
	    // return _test[_i][_qp]*(_coupled_grad_x[_qp](0) * _u[_qp] + _grad_u[_qp](0) * _coupled_x[_qp]);
	   return -_grad_test[_i][_qp](0)*_coupled_x[_qp] * _u[_qp];

	return 0.0;
}

Real
Advection_c::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian

    if(_dim == 2)
    	return -(_grad_test[_i][_qp](0)*_coupled_x[_qp] + _grad_test[_i][_qp](1)*_coupled_y[_qp])*_phi[_j][_qp];
    else
	//return _test[_i][_qp]*(_coupled_grad_x[_qp](0) * _phi[_j][_qp] + _grad_phi[_j][_qp](0) * _coupled_x[_qp]);
       return -_grad_test[_i][_qp](0)*_coupled_x[_qp]*_phi[_j][_qp];


    return 0.0;
}

Real
Advection_c::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
    if(_dim == 2){
    	int _coupled_x_aux = coupled("coupled_x");//check if this is the right component
    	int _coupled_y_aux = coupled("coupled_y");//check if this is the right component
    	if (jvar==_coupled_x_aux)   return -(_grad_test[_i][_qp](0)*_phi[_j][_qp])*_u[_qp];
    	if (jvar==_coupled_y_aux)   return -(_grad_test[_i][_qp](1)*_phi[_j][_qp])*_u[_qp];
    }else{
	int _coupled_x_aux = coupled("coupled_x");//check if this is the right component
	if (jvar==_coupled_x_aux)   
	//return _test[_i][_qp]*(_grad_phi[_j][_qp](0) * _u[_qp] + _grad_u[_qp](0) * _phi[_j][_qp]);
		
	return -(_grad_test[_i][_qp](0)*_phi[_j][_qp])*_u[_qp];

    }

	return 0.0;
}






