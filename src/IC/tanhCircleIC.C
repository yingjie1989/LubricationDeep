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

#include "tanhCircleIC.h"


template<>
InputParameters validParams<tanhCircleIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<int>("dim","dimension");
  params.addRequiredParam<std::vector<Real> >("pos", "position of the source");
  params.addRequiredParam<Real>("m_param","mparam");
  params.addRequiredParam<Real>("A_param","Aparam");
  params.addParam<Real>("outvalue",0.0,"outvalue");
  params.addParam<Real>("invalue",1.0,"invalue");
  params.addParam<Real>("radius",1.0,"invalue");

 
  return params;
}


tanhCircleIC::tanhCircleIC(const InputParameters & parameters) :
   InitialCondition(parameters),
     _dim(getParam<int>("dim")),
     _pos(getParam<std::vector<Real> >("pos")),
     _invalue(getParam<Real>("invalue")),
     _outvalue(getParam<Real>("outvalue")),
     _rad(getParam<Real>("radius")),
     _m(getParam<Real>("m_param")),
     _A(getParam<Real>("A_param"))
{}

Real
tanhCircleIC::value(const Point & p)
{
	Real r0;

	if(_dim == 1)
		r0 = std::abs(p(0) - _pos[0]);
	else
		r0 = std::sqrt( (p(0)-_pos[0]) * (p(0)-_pos[0]) + (p(1)-_pos[1]) * (p(1)-_pos[1]) );

	return _m/2*( 1 - std::tanh(_A*(r0-_rad)) ) * ( _invalue - _outvalue ) + _outvalue;

}

