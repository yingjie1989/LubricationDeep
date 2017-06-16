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

#include "SourceGradient.h"


template<>
InputParameters validParams<SourceGradient>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled",0, "Coupled variable for component x");
  params.addCoupledVar("coupled_old",0, "Coupled variable for component x"); 
  params.addCoupledVar("coupled_grad",0, "Coupled variable for component x");
  params.addParam<Real>("coeff",1.0, "coefficient");


  return params;
}


SourceGradient::SourceGradient(const InputParameters & parameters) :
   Kernel(parameters),
    _coeff(getParam<Real>("coeff")),
    _coupled(coupledValue("coupled")),
    _coupled_old(coupledValueOld("coupled_old")),
    _coupled_grad(coupledGradient("coupled_grad"))
{}

Real
SourceGradient::computeQpResidual()
{
    return -_grad_test[_i][_qp] * _coeff * _coupled_grad[_qp];
}

Real
SourceGradient::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian

	return 0.0;

}

Real
SourceGradient::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}






