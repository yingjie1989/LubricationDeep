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

#include "InCompModify.h"


template<>
InputParameters validParams<InCompModify>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled_u",0, "Coupled velocity component");   
  params.addCoupledVar("coupled_v",0, "Coupled velocity component");  
  params.addParam<int>("direct", "dimensional component");
  params.addParam<Real>("viscosity",0.0, "the viscous regularity");
 
  return params;
}


InCompModify::InCompModify(const InputParameters & parameters) :
   Kernel(parameters),
   _dim (getParam<int>("direct")),
   _mu(getParam<Real>("viscosity")),
   _coupled_gradv(coupledGradient("coupled_v")),
   _coupled_gradu(coupledGradient("coupled_u"))

{}

Real
InCompModify::computeQpResidual()
{
      return _test[_i][_qp] * (_coupled_gradu[_qp](0) + _coupled_gradv[_qp](1) + _grad_u[_qp](_dim) + _mu * _u[_qp] );

}

Real
InCompModify::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian
             
      return _test[_i][_qp] * ( _grad_phi[_j][_qp](_dim) + _mu * _phi[_j][_qp] );
}

Real
InCompModify::computeQpOffDiagJacobian(unsigned int jvar)
{
 
  // get the coupled variable jvar is referring to
        int _coupled_aux_u = coupled("coupled_u");//check if this is the right component
    	int _coupled_aux_v = coupled("coupled_v");//check if this is the right component
	if( jvar == _coupled_aux_u )
		return _test[_i][_qp] * _grad_phi[_j][_qp](0);
	else if( jvar == _coupled_aux_v )
	        return _test[_i][_qp] * _grad_phi[_j][_qp](1);
        else
		return 0.0;
}


