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

#include "ModifyReaction.h"

template<>
InputParameters validParams<ModifyReaction>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("timeStep",0,"timeStep");
  return params;
}

ModifyReaction::ModifyReaction(const InputParameters & parameters) :
    Kernel(parameters),
    _timeStep(getParam<Real>("timeStep"))

{}

Real
ModifyReaction::computeQpResidual()
{
  return 1.0/_timeStep * _test[_i][_qp] * _u[_qp];
}

Real
ModifyReaction::computeQpJacobian()
{
  return 1.0/_timeStep * _test[_i][_qp] * _phi[_j][_qp];
}

