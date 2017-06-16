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

#include "SurfPresAux.h"
#include "MooseMesh.h"
#include "libmesh/fe_interface.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"


template<>
InputParameters validParams<SurfPresAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("viscosity", 0.0, "Scalar value used for our auxiliary calculation");
  params.addParam<Real>("gamma0", 1.0, "Scalar value used for our auxiliary calculation");
  params.addParam<unsigned>("dim", 1, "dimension");
  params.addCoupledVar("u", 0.0, "Coupled variable u ");
  params.addCoupledVar("v", 0.0, "Coupled variable v ");
  params.addCoupledVar("w", 0.0, "Coupled variable w ");
  params.addCoupledVar("coupled_c",0.0, "Flux Value at the Boundary");
  params.addCoupledVar("coupled_h",0.0, "Flux Value at the Boundary");

  return params;
}

SurfPresAux::SurfPresAux(const InputParameters & parameters) :
AuxKernel(parameters),

  _coupled_c(coupledValue("coupled_c")),
//  _h_number(coupled("coupled_h")),
  _grad_h(coupledGradient("coupled_h")),
  _h_second(coupledSecond("coupled_h")),
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),
  _dim(getParam<unsigned>("dim")),
  _gamma0(getParam<Real>("gamma0")),
  _mu(getParam<Real>("viscosity"))

  
	
{}


Real
SurfPresAux::computeValue()
{

	Real norm2d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) );
	RealVectorValue normal2d(-_grad_h[_qp](0)/norm2d, 1/norm2d, 0);
	
	Real norm3d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );

	RealVectorValue normal3d(-_grad_h[_qp](0)/norm3d, -_grad_h[_qp](1)/norm3d, 1/norm3d);

	Real DUnDn = 0.0;
	Real gamma; Real Rinv;
	if(_dim == 1){
		Real dudn = normal2d * _grad_u_vel[_qp]; // + normal2d(1) * _grad_u_vel(1);  
		Real dvdn = normal2d * _grad_v_vel[_qp];	
		RealVectorValue DUDn(dudn, dvdn, 0.0);
		DUnDn = DUDn * normal2d;
        	gamma = _gamma0 * ( 1 - _coupled_c[_qp] );

				
		
		
        	//Rinv = _h_second_dev/(norm2d*norm2d*norm2d);
	}else{
		Real dudn = normal3d * _grad_u_vel[_qp];  
		Real dvdn = normal3d * _grad_v_vel[_qp]; 
		Real dwdn = normal3d * _grad_w_vel[_qp];
		RealVectorValue DUDn(dudn, dvdn, dwdn);
		DUnDn = DUDn * normal3d;
		gamma = _gamma0 * ( 1 - _coupled_c[_qp] );
                Rinv = ( _h_second[_qp](0,0) + _h_second[_qp](1,1) )/(norm3d*norm3d*norm3d);	
	}
	

	return 2.0 * _mu * DUnDn - gamma * Rinv; 

}
