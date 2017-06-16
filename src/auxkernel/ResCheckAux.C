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

#include "ResCheckAux.h"

template<>
InputParameters validParams<ResCheckAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("value", 0.0, "Scalar value used for our auxiliary calculation");
  params.addCoupledVar("coupled_u", 0.0, "Coupled variable u ");
  params.addCoupledVar("coupled_v", 0.0, "Coupled variable v");
  params.addCoupledVar("coupled_w", 0.0, "Coupled variable w");
  return params;
}

ResCheckAux::ResCheckAux(const InputParameters & parameters) :
AuxKernel(parameters),

 _coupled_u(coupledValue("coupled_u")),
 _coupled_v(coupledValue("coupled_v")),
 _coupled_w(coupledValue("coupled_w")),
 _coupled_gradu(coupledGradient("coupled_u")),
 _coupled_gradv(coupledGradient("coupled_v")),
 _coupled_gradw(coupledGradient("coupled_w")),
  _value(getParam<Real>("value"))

{}


Real
ResCheckAux::computeValue()
{
  return _coupled_gradu[_qp](0) + _coupled_gradv[_qp](1) + _coupled_gradw[_qp](2) - _value;
}
