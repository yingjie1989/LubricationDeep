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

#include "SourceMaterial.h"

//    function = '1e4/(1 + 0.01*cos(50*pi*x)*cos(50*pi*y) )'


template<>
InputParameters validParams<SourceMaterial>()
{
  InputParameters params = validParams<Material>();

  // Vectors for Linear Interpolation
  params.addRequiredParam<Real>("coeff", "Coefficient of perturbation");

  return params;
}

SourceMaterial::SourceMaterial(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _diffusivity(declareProperty<Real>("diffusivity")),
     _coeff(getParam<Real>("coeff"))

{}

void
SourceMaterial::computeQpProperties()
{

   Real u_N = _coeff*( _q_point[_qp](0) * _q_point[_qp](0) + _q_point[_qp](1) * _q_point[_qp](1) );

   Real w_N = -_coeff*( 2*_q_point[_qp](0)*_q_point[_qp](1) + _q_point[_qp](0)*_q_point[_qp](0) );

  _diffusivity[_qp] = u_N/0.01 + u_N * _coeff * 2 * _q_point[_qp](0) + w_N * _coeff * 2 * _q_point[_qp](1) - 0.01 * _coeff * 2;

}
