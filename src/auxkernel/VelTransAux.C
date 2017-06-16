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

#include "VelTransAux.h"

template<>
InputParameters validParams<VelTransAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("value", 0.0, "Scalar value used for our auxiliary calculation");
  params.addCoupledVar("coupled", 0.0, "Coupled variable u ");
  return params;
}

VelTransAux::VelTransAux(const InputParameters & parameters) :
AuxKernel(parameters),

 _coupled(coupledValue("coupled")),
  _value(getParam<Real>("value"))

{}


Real
VelTransAux::computeValue()
{
  return  _coupled[_qp];
}
