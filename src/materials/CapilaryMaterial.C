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

#include "CapilaryMaterial.h"

//    function = '1e4/(1 + 0.01*cos(50*pi*x)*cos(50*pi*y) )'


template<>
InputParameters validParams<CapilaryMaterial>()
{
  InputParameters params = validParams<Material>();

  // Vectors for Linear Interpolation
  params.addRequiredParam<Real>("kappa0", "The constant capilary coefficient");
  params.addRequiredParam<Real>("eta0", "The perturbation coefficient");
  params.addRequiredParam<std::vector<Real> >("coeff", "Coefficient of perturbation");
  params.addRequiredParam<std::vector<Real> >("wavex", "The wavenumber in x direction");
  params.addRequiredParam<std::vector<Real> >("wavey", "The wavenumber in y direction");

//  params.addCoupledVar("diffusion_gradient", "The gradient of this variable will be used to compute a velocity vector property.");

  return params;
}

CapilaryMaterial::CapilaryMaterial(const InputParameters & parameters) :
    Material(parameters),
    // Declare that this material is going to provide a Real
    // valued property named "diffusivity" that Kernels can use.
    _diffusivity(declareProperty<Real>("diffusivity")),

     _kappa0(getParam<Real>("kappa0")),
     _eta0(getParam<Real>("eta0")),
     _coeff(getParam<std::vector<Real> >("coeff")),
     _wavex(getParam<std::vector<Real> >("wavex")),
     _wavey(getParam<std::vector<Real> >("wavey"))

{}

void
CapilaryMaterial::computeQpProperties()
{
 
  Real sum = 0.0;
  for(unsigned int i=0; i<_wavex.size(); i++) 
	sum += _coeff[i]*(std::cos(_wavex[i]*_q_point[_qp](0)) * std::cos(_wavey[i]*_q_point[_qp](1)) + 1);

  _diffusivity[_qp] = _kappa0 + _eta0 * sum * std::exp(-5*(_q_point[_qp](0)-1)*(_q_point[_qp](0)-1) );

}
