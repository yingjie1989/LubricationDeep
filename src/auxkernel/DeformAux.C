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

#include "DeformAux.h"

template<>
InputParameters validParams<DeformAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("value", 0.0, "Scalar value used for our auxiliary calculation");
  params.addRequiredCoupledVar("coupled", "Coupled variable");
  return params;
}

DeformAux::DeformAux(const InputParameters & parameters) :
AuxKernel(parameters),

 _coupled_val(coupledValue("coupled")),

  _value(getParam<Real>("value"))

{}


Real
DeformAux::computeValue()
{
  return _coupled_val[_qp] * ( _q_point[_qp](1) - _value )/std::abs(_value);
}
