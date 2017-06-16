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

#include "Advection_h.h"


template<>
InputParameters validParams<Advection_h>()
{
  InputParameters params = validParams<Diffusion>();
  params.addCoupledVar("coupled_x", "Coupled variable for component x");
  params.addCoupledVar("coupled_y", "Coupled variable for component y");
  //params.addRequiredParam<std::string>("mat_name", "Praft material name");
  params.addParam<MaterialPropertyName>("f_param","fmat","");
  return params;
}


Advection_h::Advection_h(const InputParameters & parameters) :
   Diffusion(parameters),
    //_mat_name(getParam<std::string>("mat_name")),
    _fmat(getMaterialProperty<Real>("f_param")),
    _coupled_x(coupledValue("coupled_x")),
    _coupled_y(coupledValue("coupled_y"))

{}

Real
Advection_h::computeQpResidual()
{

        return -(_grad_test[_i][_qp](0)*_coupled_x[_qp] + _grad_test[_i][_qp](1)*_coupled_y[_qp])*(_u[_qp] - _fmat[_qp]);
}

Real
Advection_h::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian

    return -(_grad_test[_i][_qp](0)*_coupled_x[_qp] + _grad_test[_i][_qp](1)*_coupled_y[_qp])*_phi[_j][_qp];


}

Real
Advection_h::computeQpOffDiagJacobian(unsigned int jvar)
{
 
  // get the coupled variable jvar is referring to
    int _coupled_x_aux = coupled("coupled_x");//check if this is the right component
    int _coupled_y_aux = coupled("coupled_y");//check if this is the right component
    if (jvar==_coupled_x_aux)   return -(_grad_test[_i][_qp](0)*_phi[_j][_qp])*(_u[_qp] - _fmat[_qp]);
    if (jvar==_coupled_y_aux)   return -(_grad_test[_i][_qp](1)*_phi[_j][_qp])*(_u[_qp] - _fmat[_qp]);

	return 0.0;
}






