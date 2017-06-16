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

#include "FingerCircleIC.h"


template<>
InputParameters validParams<FingerCircleIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("b_param","bparam");
  params.addRequiredParam<Real>("A_param","Aparam");
  params.addRequiredParam<Real>("B_param","Bparam");
  params.addParam<Real>("C_param","Cparam");
  params.addParam<std::vector<Real> >("coeff","coefficient for underlying fluid");
  params.addRequiredParam<bool>("hIC","varIndicator");
  return params;
}


FingerCircleIC::FingerCircleIC(const InputParameters & parameters) :
   InitialCondition(parameters),
     _b(getParam<Real>("b_param")),
     _A(getParam<Real>("A_param")),
     _B(getParam<Real>("B_param")),
     _C(getParam<Real>("C_param")), 
    _coeff(getParam<std::vector<Real> >("coeff")),
    _hIC(getParam<bool>("hIC"))
{}

Real
FingerCircleIC::value(const Point & p)
{
	
	if(_hIC){
		Real hdrop = 0.0;
	        Real r0 = std::sqrt( p(0)*p(0) + p(1)*p(1) );
		Real sum = 0.0;
		for(unsigned int i=0; i<_coeff.size(); i++)
			sum += _coeff[i];
		 
		hdrop =  ( 1 - r0*r0 + _b)*funcH(1-r0) + _b*funcH( r0-1 ) + sum * _C; // * std::exp(-_B*( r0 - 1 )*( r0 - 1 ) );

	  	Real hpert = 0.0;

		return hdrop + (-_A + 2*_A*(double)rand()/RAND_MAX )*std::exp(-_B*( r0 - 1 )*( r0 - 1 ) );
	}else{
		Real r0 = std::sqrt( p(0)*p(0) + p(1)*p(1) );
		return funcH( 1 - r0 );
	}

	return 0.0;

}

Real 
FingerCircleIC::funcH(Real x)
{
	
	return 0.5*( 1 + ( std::exp(20*x) - std::exp(-20*x) )/( std::exp(20*x) + std::exp(-20*x) ) );
}
