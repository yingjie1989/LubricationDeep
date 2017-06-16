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

#include "SourceManuFac.h"


template<>
InputParameters validParams<SourceManuFac>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<MaterialPropertyName>("f_param","fmat","");
  

  return params;
}


SourceManuFac::SourceManuFac(const InputParameters & parameters) :
   Kernel(parameters),
   _fmat(getMaterialProperty<Real>("f_param"))
   
{}

Real
SourceManuFac::computeQpResidual()
{
    return -_test[_i][_qp] * _fmat[_qp];
}

Real
SourceManuFac::computeQpJacobian()
{
  //remember that this part computes the diagonal part of jacobian

	return 0.0;

}

Real
SourceManuFac::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}






