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

#include "ReactionCoupled.h"

template<>
InputParameters validParams<ReactionCoupled>()
{
  InputParameters params = validParams<Kernel>();

  params.addParam<MaterialPropertyName>("rho_param", "rho", "density");
  params.addParam<MaterialPropertyName>("grav_param", "grav", "gravity");
  params.addCoupledVar("coupled", "Coupled variable, height");
  params.addParam<Real>("coeff","coefficient");

  
  return params;
}

ReactionCoupled::ReactionCoupled(const InputParameters & parameters) :
    Kernel(parameters),
     _rho(getMaterialProperty<Real>("rho_param")),
    _grav(getMaterialProperty<Real>("grav_param")),
    _coeff(getParam<Real>("coeff")),
    _coupled_h(coupledValue("coupled"))
{}

Real
ReactionCoupled::computeQpResidual()
{

	 return -_coeff*_rho[_qp]*_grav[_qp]*_test[_i][_qp] * _coupled_h[_qp];

}

Real
ReactionCoupled::computeQpJacobian()
{
	return 0.0;
}

Real
ReactionCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
    int _coupled_h_aux = coupled("coupled");
    if (jvar==_coupled_h_aux)
        return -_coeff*_rho[_qp]*_grav[_qp]*_test[_i][_qp]*_phi[_j][_qp];

    return 0.0;

}


