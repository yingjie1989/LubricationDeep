/****************************************************************/
/*               DO NOT MODIFY THIS CLASS                       */
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

#include "NSPressure.h"

template<>
InputParameters validParams<NSPressure>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<MaterialPropertyName>("grav_param", "grav", "gravity");
  params.addCoupledVar("coupled", "Coupled variable, height");
  params.addParam<int>("direct","coefficient");
  params.addParam<Real>("coeff",1.0,"coefficient");
  
  return params;
}

NSPressure::NSPressure(const InputParameters & parameters) :
    Kernel(parameters),
    _grav(getMaterialProperty<Real>("grav_param")),
    _dim(getParam<int>("direct")),
    _coeff(getParam<Real>("coeff")),
    _coupled(coupledValue("coupled")),
    _coupled_grad(coupledGradient("coupled"))
{}

Real
NSPressure::computeQpResidual()
{
	 return _test[_i][_qp] * _grav[_qp] - _coeff * _grad_test[_i][_qp](_dim) * _coupled[_qp];
}

Real
NSPressure::computeQpJacobian()
{
	return 0.0;
}

Real
NSPressure::computeQpOffDiagJacobian(unsigned int jvar)
{
    int _coupled_aux = coupled("coupled");
    if (jvar==_coupled_aux)
	return -_coeff * _grad_test[_i][_qp](_dim) * _phi[_j][_qp];

    return 0.0;

}


