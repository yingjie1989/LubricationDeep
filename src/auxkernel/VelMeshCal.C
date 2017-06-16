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

#include "VelMeshCal.h"

template<>
InputParameters validParams<VelMeshCal>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("timeStep", 1.0, "Scalar value used for our auxiliary calculation");
  params.addCoupledVar("coupled", 0.0, "Coupled variable u ");
  return params;
}

VelMeshCal::VelMeshCal(const InputParameters & parameters) :
AuxKernel(parameters),

 _coupled(coupledValue("coupled")),
 _coupledOld(coupledValueOld("coupled")),
 _DeltaT(getParam<Real>("timeStep"))

{}


Real
VelMeshCal::computeValue()
{
  return  ( _coupled[_qp] - _coupledOld[_qp] )/_DeltaT;
}
