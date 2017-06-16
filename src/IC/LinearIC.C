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

#include "LinearIC.h"

template<>
InputParameters validParams<LinearIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<int>("direct", "The direction of variation");
  params.addRequiredParam<Real>("coefficient", "The value of the initial condition");
      return params;
}

LinearIC::LinearIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _dim(getParam<int>("direct")),
    _coefficient(getParam<Real>("coefficient"))
	
{}

Real
LinearIC::value(const Point & p)
{
    return _coefficient*std::abs(p(0))*std::abs(p(1));
}
