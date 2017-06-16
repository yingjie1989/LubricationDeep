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

#include "InCompressVelocity.h"


template<>
InputParameters validParams<InCompressVelocity>()
{
  InputParameters params = validParams<Diffusion>();
//  params.addRequiredParam<std::string>("mat_name", "material property");
  params.addCoupledVar("coupled", "Coupled velocity component");
  params.addParam<int>("direct", "dimensional component");

  return params;
}


InCompressVelocity::InCompressVelocity(const InputParameters & parameters) :
   Diffusion(parameters),
   
   _dim (getParam<int>("direct")),
  // _coupled(coupledValue("coupled"))
   _coupled_gradv(coupledGradient("coupled"))

{}

Real
InCompressVelocity::computeQpResidual()
{

	return _test[_i][_qp] * ( _u[_qp] + _coupled_gradv[_qp](_dim) );
}

Real
InCompressVelocity::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian
	 return _test[_i][_qp] * ( _phi[_j][_qp] );

}

Real
InCompressVelocity::computeQpOffDiagJacobian(unsigned int jvar)
{
 
  // get the coupled variable jvar is referring to

    	int _coupled_aux = coupled("coupled");//check if this is the right component
	if( jvar == _coupled_aux )
		return _test[_i][_qp] * _grad_phi[_j][_qp](_dim);


	return 0.0;
}






