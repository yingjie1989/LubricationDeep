/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeSurfaceStress.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<ComputeSurfaceStress>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<AuxVariableName>("surface_disp", "The enrichment displacement");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  params.addParam<std::vector<Real> >("Reyold", "Value multiplied by the coupled value on the boundary");
  params.addParam<Real>("mesh_size", "indicate tolerance");
  params.addParam<Real>("tol", "indicate tolerance");
  params.addParam<unsigned>("direct", "indicate dimension");
  params.addRequiredCoupledVar("coupled_c", "Flux Value at the Boundary");
  params.addCoupledVar("coupled_h", "Flux Value at the Boundary");
  
return params;
}

ComputeSurfaceStress::ComputeSurfaceStress(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _surface_stress(declareProperty<RealVectorValue>("surface_stress")),
    _hele(getParam<Real>("mesh_size")),
    _Re(getParam<std::vector<Real> >("Reyold")),
    _tol(getParam<Real>("tol")),
    _dim (getParam<unsigned>("direct")),
    _some_var_val(coupledValue("coupled_c")),
    _c_grad(coupledGradient("coupled_c")),
    _grad_h(coupledGradient("coupled_h")),
   // _h_var_number(coupled("coupled_h")),
    _h_second(coupledSecond("coupled_h"))
{

    //const AuxVariableName & au_vnames = getParam<AuxVariableName>("surface_disp");
    //AuxiliarySystem & au = _fe_problem.getAuxiliarySystem();
    //_h_var_number = &(au.getVariable(0, au_vnames));   

}

void
ComputeSurfaceStress::computeQpProperties()
{
  /*FEType fe_type(Utility::string_to_enum<Order>("first"),Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule (_qrule);
  
  fe->reinit (_current_elem); 
  AuxiliarySystem & au = _fe_problem.getAuxiliarySystem();
  const NumericVector<Number> & auxsln = *au.currentSolution();*/
 
        
	Real norm2d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) );
   	RealVectorValue normal2d(-_grad_h[_qp](0)/norm2d, 1/norm2d, 0);
   	RealVectorValue tau2d(1/norm2d,_grad_h[_qp](0)/norm2d,0);

	Real norm3d = sqrt( 1 + _grad_h[_qp](0) * _grad_h[_qp](0) + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue normal3d(-_grad_h[_qp](0)/norm3d, -_grad_h[_qp](1)/norm3d, 1/norm3d);
   Real norm2dy = sqrt( 1 + _grad_h[_qp](1) * _grad_h[_qp](1) );
   RealVectorValue tau3dx(1/norm2d, 0, _grad_h[_qp](0)/norm2d);
   RealVectorValue tau3dy(0, 1/norm2dy, _grad_h[_qp](1)/norm2dy);

   	Real sigma_val = 0.0;
   	Real _sigma_n = 0.0;
   	Real _sigma_tx = 0.0;
   	Real _sigma_ty = 0.0;


	/*for(unsigned i=0; i<_current_elem->n_nodes(); i++)
	{
		Node * node_i = _current_elem->get_node(i);
        dof_id_type dofi = node_i->dof_number(au.number(),_h_var_number->number(), 0);
        Real soln_locali = auxsln(dofi);
		printf("solution %lf, ",soln_locali);
	}*/

	if(_some_var_val[_qp] > _tol){
	
	   if(_dim == 1){
        _sigma_tx = -_Re[0]*( _c_grad[_qp](0) * tau2d(0) + _c_grad[_qp](1) * tau2d(1) );

        Real gamma = _Re[1]*( 1 - _some_var_val[_qp] );
   
	/*Node * node_0 = _current_elem->get_node(3);
        dof_id_type dof0 = node_0->dof_number(au.number(),_h_var_number->number(), 0);
        Real soln_local0 = auxsln(dof0);
   
        Node * node_1 = _current_elem->get_node(6);
        dof_id_type dof1 = node_1->dof_number(au.number(),_h_var_number->number(), 0);
	Real soln_local1 = auxsln(dof1);

        Node * node_2 = _current_elem->get_node(2);
        dof_id_type dof2 = node_2->dof_number(au.number(),_h_var_number->number(), 0);
	Real soln_local2 = auxsln(dof2);

        Real _h_second_00 = ( soln_local0 + soln_local2 - 2*soln_local1 )/_hele/_hele;

        Real Rinv = _h_second_00/(norm2d*norm2d*norm2d);*/

	Real Rinv = _h_second[_qp](0,0)/(norm2d*norm2d*norm2d);

        _sigma_n  = gamma * Rinv;
	
   }else if(_dim == 2){

	_sigma_tx = -_Re[0]*( _c_grad[_qp] * tau3dx );
          _sigma_ty = -_Re[1]*( _c_grad[_qp] * tau3dy );

        Real gamma = _Re[2] * (1 - _some_var_val[_qp]);
        Real Rinv = 0.5*( (1 + _grad_h[_qp](0)*_grad_h[_qp](0))*_h_second[_qp](1,1) + (1 + _grad_h[_qp](1)*_grad_h[_qp](1))*_h_second[_qp](0,0) - 2*_grad_h[_qp](0)*_grad_h[_qp](1)*_h_second[_qp](0,1) )/(norm3d*norm3d*norm3d);

        _sigma_n  = gamma * Rinv;
    }}

	RealVectorValue surfaceT(_sigma_tx,_sigma_ty,_sigma_n);

	_surface_stress[_qp] = surfaceT;
}

void
ComputeSurfaceStress::initQpStatefulProperties()
{
}
