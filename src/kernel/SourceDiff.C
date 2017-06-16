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

#include "SourceDiff.h"


template<>
InputParameters validParams<SourceDiff>()
{
  InputParameters params = validParams<Diffusion>();
  params.addCoupledVar("coupled",0, "Coupled variable for component x");
  params.addCoupledVar("coupled_old",0, "Coupled variable for component x"); 
  params.addCoupledVar("coupled_grad",0, "Coupled variable for component x");
  params.addParam<int>("component",0, "gradient component");
  params.addParam<Real>("coeff",1.0, "coefficient");


  return params;
}


SourceDiff::SourceDiff(const InputParameters & parameters) :
   Diffusion(parameters),
    _component(getParam<int>("component")),
    _coeff(getParam<Real>("coeff")),
    _coupled(coupledValue("coupled")),
    _coupled_old(coupledValueOld("coupled_old")),
    _coupled_grad(coupledGradient("coupled_grad"))
{}

Real
SourceDiff::computeQpResidual()
{
    return -_test[_i][_qp] * _coeff * ( _coupled[_qp] + _coupled_grad[_qp](_component) + _coupled_old[_qp] );
}

Real
SourceDiff::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian

	return 0.0;

}

Real
SourceDiff::computeQpOffDiagJacobian(unsigned int jvar)
{
	int _coupled_aux = coupled("coupled");
	int _coupled_grad_aux = coupled("coupled_grad");

	if(jvar == _coupled_aux)
		return -_test[_i][_qp] * _coeff * _phi[_j][_qp];
	//if(jvar == _coupled_grad_aux)
        //        return -_test[_i][_qp] * _grad_phi[_j][_qp](_component);
	return 0.0;
}






