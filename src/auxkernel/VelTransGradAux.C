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

#include "VelTransGradAux.h"

template<>
InputParameters validParams<VelTransGradAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("value", 0.0, "Scalar value used for our auxiliary calculation");
  params.addParam<int>("component", 0, "Scalar value used for our auxiliary calculation");
  params.addCoupledVar("phi", 0.0, "Coupled variable u ");
  params.addCoupledVar("coupled", 0.0, "Coupled variable u ");
  return params;
}

VelTransGradAux::VelTransGradAux(const InputParameters & parameters) :
AuxKernel(parameters),

 _coupled_grad(coupledGradient("coupled")),
  _field(coupledValue("phi")),
  _component(getParam<int>("component")),
  _value(getParam<Real>("value"))

{}


Real
VelTransGradAux::computeValue()
{
  return  _value * _coupled_grad[_qp](_component) * _field[_qp];
}
