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

#include "CoupledNeumannStressBC.h"

template<>
InputParameters validParams<CoupledNeumannStressBC>()
{
	InputParameters params = validParams<IntegratedBC>();
          params.addParam<MaterialPropertyName>("surface_stress","sursigma","");
	params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<unsigned>("direct", "indicate dimension");
        params.addParam<unsigned>("component", "indicate component");
        params.addRequiredCoupledVar("coupled_c", "Flux Value at the Boundary");
        params.addRequiredCoupledVar("coupled_h", "Flux Value at the Boundary");

        return params;
}

CoupledNeumannStressBC::CoupledNeumannStressBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _surface_stress(getMaterialProperty<RealVectorValue>("surface_stress")),
    _dim (getParam<unsigned>("direct")),
    _component(getParam<unsigned>("component")), 
    _some_var_val(coupledValue("coupled_c")),
    _c_grad(coupledGradient("coupled_c")),
    _grad_h(coupledGradient("coupled_h")) 
     
{}

Real
CoupledNeumannStressBC::computeQpResidual()
{

   Real norm2d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) );
   RealVectorValue normal2d(-_grad_h[_qp](0)/norm2d, 1/norm2d, 0);
   RealVectorValue tau2d(1/norm2d,_grad_h[_qp](0)/norm2d,0);

   Real norm3d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue normal3d(-_grad_h[_qp](0)/norm3d, -_grad_h[_qp](1)/norm3d, 1/norm3d);
   Real normp = sqrt( _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue tau3d(_grad_h[_qp](0)/norm3d/normp, _grad_h[_qp](1)/norm3d/normp, normp/norm3d);

   Real sigma_val = 0.0;
   Real _sigma_n = _surface_stress[_qp](2);
   Real _sigma_tx = _surface_stress[_qp](0);
   Real _sigma_ty = _surface_stress[_qp](1);


   if(_dim == 1){

	if(_component == 0 )
		sigma_val = _sigma_tx; // * tau2d(0) + _sigma_n * normal2d(0);
	else if(_component == 1)
		sigma_val = _sigma_n; // * normal2d(1) + _sigma_tx * tau2d(1);
	else
		sigma_val = 0.0;

   }else if(_dim == 2){ 
        
 
	if(_component == 0)
		sigma_val = _sigma_tx; 
	else if(_component == 1)
		sigma_val = _sigma_ty;
	else if(_component == 2)
		sigma_val = _sigma_n;
	else
		sigma_val =  0.0;
   }
	return - _test[_i][_qp] * sigma_val;
}
