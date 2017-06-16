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

#include "QuadraticIC.h"

template<>
InputParameters validParams<QuadraticIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<int>("direct", "The direction of variation");
  params.addRequiredParam<Real>("coefficient", "The value of the initial condition");
      return params;
}

QuadraticIC::QuadraticIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _dim(getParam<int>("direct")),
    _coefficient(getParam<Real>("coefficient"))
	
{}

Real
QuadraticIC::value(const Point & p)
{
    return _coefficient*( 25 - (p(_dim) - 5.0)*(p(_dim) - 5.0) );
}
