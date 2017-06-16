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

#include "SourceMaterial3D.h"

//    function = '1e4/(1 + 0.01*cos(50*pi*x)*cos(50*pi*y) )'


template<>
InputParameters validParams<SourceMaterial3D>()
{
  InputParameters params = validParams<Material>();

  // Vectors for Linear Interpolation
  params.addRequiredParam<Real>("coeff", "Coefficient of perturbation");
  params.addParam<Real>("timestep", 0.01, "timestep");
  params.addParam<Real>("mu", 0.01, "timestep");
  params.addParam<Real>("direct", 0, "component");

  return params;
}

SourceMaterial3D::SourceMaterial3D(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _diffx(declareProperty<Real>("diffusex")),
    _diffy(declareProperty<Real>("diffusey")),
     _coeff(getParam<Real>("coeff")),
    _dt(getParam<Real>("timestep")),
    _mu(getParam<Real>("mu")),
    _component(getParam<Real>("direct"))

{}

void
SourceMaterial3D::computeQpProperties()
{

   Real u_N = _coeff*( _q_point[_qp](0) * _q_point[_qp](0)  + _q_point[_qp](1) * _q_point[_qp](1) +  _q_point[_qp](2) * _q_point[_qp](2) );

   Real v_N = _coeff*( _q_point[_qp](0) * _q_point[_qp](0) + _q_point[_qp](1) * _q_point[_qp](1)  +  _q_point[_qp](2) * _q_point[_qp](2) );


   Real w_N = -_coeff*( 2*_q_point[_qp](0)*_q_point[_qp](2) + 2*_q_point[_qp](1)*_q_point[_qp](2) + _q_point[_qp](0)*_q_point[_qp](0) + _q_point[_qp](1)*_q_point[_qp](1) );


   
        _diffx[_qp] = u_N/_dt + u_N * _coeff * 2 * _q_point[_qp](0) + v_N * _coeff * 2 * _q_point[_qp](1) + w_N * _coeff * 2 * _q_point[_qp](2) - _mu * _coeff * 2;
  

        _diffy[_qp] = v_N/_dt + u_N * _coeff * 2 * _q_point[_qp](0) + v_N * _coeff * 2 * _q_point[_qp](1) + w_N * _coeff * 2 * _q_point[_qp](2) - _mu * _coeff * 2;

}
