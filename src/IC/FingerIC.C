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

#include "FingerIC.h"


template<>
InputParameters validParams<FingerIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("b_param","bparam");
  params.addRequiredParam<Real>("A_param","Aparam");
  params.addRequiredParam<Real>("B_param","Bparam");
  params.addRequiredParam<bool>("hIC","varIndicator");
  return params;
}


FingerIC::FingerIC(const InputParameters & parameters) :
   InitialCondition(parameters),
     _b(getParam<Real>("b_param")),
     _A(getParam<Real>("A_param")),
     _B(getParam<Real>("B_param")),
     _hIC(getParam<bool>("hIC"))
{}

Real
FingerIC::value(const Point & p)
{
	if(_hIC){
		Real hdrop = 0.0;
		hdrop =  ( 1 - p(0)*p(0) + _b)*funcH(1-p(0)) + _b*funcH( p(0)-1 ) + 3.5*0.01*std::exp(-_B*( p(0) - 1 )*( p(0) - 1 ) );

	  	double Ci[4] = {1.0, 1.0, 1.0, 0.5};
	  	int ki[4] = {2,5,7,20}; 
	  	Real hpert = 0.0;

		for(int i=0; i<4; i++)
			hpert += Ci[i]*(std::cos(ki[i]*p(1))+1);

		return hdrop - hpert*_A*std::exp(-_B*( p(0) - 1 )*( p(0) - 1 ) );
	}else{
		return funcH( 1 - p(0) );
	}

	return 0.0;

}

Real 
FingerIC::funcH(Real x)
{
	
	return 0.5*( 1 + ( std::exp(20*x) - std::exp(-20*x) )/( std::exp(20*x) + std::exp(-20*x) ) );
}
