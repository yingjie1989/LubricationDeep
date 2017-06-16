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

#include "NSCoupled.h"


template<>
InputParameters validParams<NSCoupled>()
{
  InputParameters params = validParams<Diffusion>();
//  params.addRequiredParam<std::string>("mat_name", "material property");
  params.addCoupledVar("coupled_v", "Coupled variable for thickness h");
  params.addParam<int>("direct", "dimensional component");

  return params;
}


NSCoupled::NSCoupled(const InputParameters & parameters) :
   Diffusion(parameters),
   _dim (getParam<int>("direct")),
   _coupled_v(coupledValue("coupled_v"))
  // _coupled_gradv(coupledGradient("coupled"))

{}

Real
NSCoupled::computeQpResidual()
{
    //printf("w is %f\n",_coupled_v[_qp]);

    return _test[_i][_qp] * _coupled_v[_qp] * _grad_u[_qp](_dim);
}

Real
NSCoupled::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian
  
  // printf("w is %f\n",_coupled_v[_qp]);

   return _test[_i][_qp] * _coupled_v[_qp] * _grad_phi[_j][_qp](_dim);
}

Real
NSCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
 
  // get the coupled variable jvar is referring to

    int _coupled_aux = coupled("coupled_v");//check if this is the right component

    if (jvar == _coupled_aux){	
        return _test[_i][_qp] * _phi[_j][_qp] * _grad_u[_qp](_dim);
     }	
     
     return 0.0;
}






