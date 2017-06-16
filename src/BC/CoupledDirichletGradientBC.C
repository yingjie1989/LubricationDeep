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

#include "CoupledDirichletGradientBC.h"

template<>
InputParameters validParams<CoupledDirichletGradientBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addParam<Real>("alpha", 1.0, "Value multiplied by the coupled value on the boundary");
  params.addParam<int>("component", 0, "Value of the gradient");
  params.addRequiredCoupledVar("some_var", "Value on the Boundary");
  return params;
}


CoupledDirichletGradientBC::CoupledDirichletGradientBC(const InputParameters & parameters) :
    NodalBC(parameters),
    _alpha(getParam<Real>("alpha")),
    _component(getParam<int>("component")),
    _some_var_grad(coupledGradient("some_var"))

{}

Real
CoupledDirichletGradientBC::computeQpResidual()
{
  return _u[_qp]-(_alpha*_some_var_grad[_qp](_component));
}
