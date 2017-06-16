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

#include "InCompressStream.h"


template<>
InputParameters validParams<InCompressStream>()
{
  InputParameters params = validParams<Diffusion>();
//  params.addRequiredParam<std::string>("mat_name", "material property");
  params.addCoupledVar("coupled", "Coupled velocity component");
  params.addParam<int>("direct", "dimensional component");

  return params;
}


InCompressStream::InCompressStream(const InputParameters & parameters) :
   Diffusion(parameters),
   
   _dim (getParam<int>("direct")),
   _coupled_v(coupledValue("coupled"))
  // _coupled_gradv(coupledGradient("coupled"))

{}

Real
InCompressStream::computeQpResidual()
{

	return _test[_i][_qp] * ( _grad_u[_qp](_dim) - _coupled_v[_qp] );
}

Real
InCompressStream::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian
	 return _test[_i][_qp] * _grad_phi[_j][_qp](_dim);

}

Real
InCompressStream::computeQpOffDiagJacobian(unsigned int jvar)
{
 
  // get the coupled variable jvar is referring to

    	int _coupled_aux = coupled("coupled");//check if this is the right component
	if( jvar == _coupled_aux )
		return -_test[_i][_qp] * _phi[_j][_qp];


	return 0.0;
}






