/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "InCompSUPG.h"

template<>
InputParameters validParams<InCompSUPG>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<MaterialPropertyName>("grav_param", "grav", "gravity");
  //Coupled Variables

  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
 
   params.addCoupledVar("u_mesh", 0, "x-meshvelocity"); // only required in 3D
  params.addCoupledVar("v_mesh", 0, "y-meshvelocity"); // only required in 3D
  params.addCoupledVar("w_mesh", 0, "z-meshvelocity"); // only required in 3D 

 // Required parameters
  params.addParam<Real>("timeStep",0.01, "pressure coefficient");
  params.addRequiredParam<std::vector<Real> >("mu", "dynamic viscosity");
  params.addParam<std::vector<Real> >("rho", "density");
  params.addParam<Real>("coeff_pressure",1.0, "pressure coefficient");
  params.addParam<Real>("coeff_param",0.0, "coefficient");
  params.addRequiredParam<std::vector<Real> >("hele", "Coefficient of perturbation");
  params.addRequiredParam<unsigned>("vert", "vertical diffusion direction");
 
  return params;
}

InCompSUPG::InCompSUPG(const InputParameters & parameters) :
    Kernel(parameters),
    _grav(getMaterialProperty<Real>("grav_param")),
   // Coupled variables
   _u_vel(coupledValue("u")),
   _v_vel(coupledValue("v")),
   _w_vel(coupledValue("w")),
  
   _u_vel_mesh(coupledValue("u_mesh")),
  _v_vel_mesh(coupledValue("v_mesh")),
  _w_vel_mesh(coupledValue("w_mesh")),

  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),
 
  // Second
  _second_u_vel(coupledSecond("u")),
  _second_v_vel(coupledSecond("v")),
  _second_w_vel(coupledSecond("w")), 
  _second_phi(secondPhi()),


  // TimeDerivative

  _u_vel_dot(coupledDot("u")),
  _v_vel_dot(coupledDot("v")),
  _w_vel_dot(coupledDot("w")),

   // TimeDerivative with respect to the coefficients

  _u_vel_dot_du(coupledDotDu("u")),
  _v_vel_dot_du(coupledDotDu("v")),
  _w_vel_dot_du(coupledDotDu("w")),


  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),

  // Required parameters
  _DeltaT(getParam<Real>("timeStep")),
  _mu(getParam<std::vector<Real> >("mu")),
  _rho(getParam<std::vector<Real> >("rho")),
  _meshsize(getParam<std::vector<Real> >("hele")),
  _coeff_pressure(getParam<Real>("coeff_pressure")),
  _coeff_param(getParam<Real>("coeff_param")),
  _vert(getParam<unsigned>("vert"))

{
}

Real
InCompSUPG::computeQpResidual()
{
  // See "Component SUPG contributions" section of notes for details.

  // Values to be summed up and returned

  // Velocity vector

   RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);


  // Velocity vector magnitude squared
  Real velmag2 = vel.norm_sq();

  Real _hele = 0.0;

  Real _muK = 0.0;

  if(_meshsize.size() == 1 ){
     _hele = _meshsize[0]; _muK = _mu[0];
  }else if(_meshsize.size() == 2){
     _hele = sqrt( _meshsize[0]*_meshsize[1] ); _muK = sqrt( _mu[0]*_mu[1] );
  }else{
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);
     _muK = std::pow(_mu[0]*_mu[1]*_mu[2],0.3333);
  }
  Real _tau_PSPG = _hele*_hele/_muK/4.0;

  Real mass_part = 0.0;

	mass_part += (_grad_test[_i][_qp](0) * _rho[0] * _u_vel_dot[_qp]
		    +_grad_test[_i][_qp](1) * _rho[1] * _v_vel_dot[_qp]
		    +_grad_test[_i][_qp](2) * _rho[2] * _w_vel_dot[_qp] );


  Real mom_part = 0.0;

     mom_part += ( _grad_test[_i][_qp](0) * _rho[0] * (vel * _grad_u_vel[_qp])
  		  +_grad_test[_i][_qp](1) * _rho[1] * (vel * _grad_v_vel[_qp])
		  +_grad_test[_i][_qp](2) * _rho[2] * (vel * _grad_w_vel[_qp]) ); 

  Real diff_part = 0.0;

  Real u_part = _second_u_vel[_qp](_vert,_vert) + _coeff_param * _second_u_vel[_qp](0,0); 
  Real v_part = _second_v_vel[_qp](_vert,_vert) + _coeff_param * _second_v_vel[_qp](0,0);
  Real w_part = _second_w_vel[_qp](_vert,_vert) + _coeff_param * _second_w_vel[_qp](0,0);

  if(_vert > 1){
    u_part += _coeff_param * _second_u_vel[_qp](1,1);        
    v_part += _coeff_param * _second_v_vel[_qp](1,1);
    w_part += _coeff_param * _second_w_vel[_qp](1,1);
  }


     diff_part -= (_mu[0] * _grad_test[_i][_qp](0) * u_part + _mu[1] *  _grad_test[_i][_qp](1) * v_part + _mu[2] * _grad_test[_i][_qp](2) * w_part );

   Real pre_part = 0.0;

   if(_vert == 1)
     pre_part += _grad_test[_i][_qp](0) * _grad_u[_qp](0)
		+_grad_test[_i][_qp](1) * ( _coeff_pressure * _grad_u[_qp](1) + _grav[_qp] );
   else
     pre_part += _grad_test[_i][_qp](0) * _grad_u[_qp](0)
                +_grad_test[_i][_qp](1) * _grad_u[_qp](1)
		+_grad_test[_i][_qp](2) * ( _coeff_pressure * _grad_u[_qp](2) + _grav[_qp] );
 
  return _tau_PSPG * ( mass_part + mom_part + diff_part + pre_part);
}

Real
InCompSUPG::computeQpJacobian()
{
	
   RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);

   
   Real velmag2 = vel.norm_sq();
   
   Real _hele = 0.0;
   Real _muK = 0.0;

  if(_meshsize.size() == 1 ){
     _hele = _meshsize[0]; _muK = _mu[0];
  }else if(_meshsize.size() == 2){
     _hele = sqrt( _meshsize[0]*_meshsize[1] ); _muK = sqrt( _mu[0]*_mu[1] );
  }else{
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);
     _muK = std::pow(_mu[0]*_mu[1]*_mu[2],0.3333);
  }
  Real _tau_PSPG = _hele*_hele/_muK/4.0;


  Real pre_part  = 0.0;

     if(_vert == 1)
	pre_part += _grad_test[_i][_qp](0) * _grad_phi[_j][_qp](0) + _coeff_pressure *  _grad_test[_i][_qp](1) * _grad_phi[_j][_qp](1);
	else
	pre_part += _grad_test[_i][_qp](0) * _grad_phi[_j][_qp](0) + _coeff_pressure *  _grad_test[_i][_qp](2) * _grad_phi[_j][_qp](2) + _grad_test[_i][_qp](1) * _grad_phi[_j][_qp](1) ;

  return _tau_PSPG * pre_part;

}

Real
InCompSUPG::computeQpOffDiagJacobian(unsigned int jvar)
{

    if (jvar == _u_vel_var_number){

	   RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);


	Real velmag2 = vel.norm_sq();

	 Real _hele = 0.0;

	Real _muK = 0.0;

  if(_meshsize.size() == 1 ){
     _hele = _meshsize[0]; _muK = _mu[0];
  }else if(_meshsize.size() == 2){
     _hele = sqrt( _meshsize[0]*_meshsize[1] ); _muK = sqrt( _mu[0]*_mu[1] );
  }else{
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);
     _muK = std::pow(_mu[0]*_mu[1]*_mu[2],0.3333);
  }
  Real _tau_PSPG = _hele*_hele/_muK/4.0;

 
     Real mass_part = 0.0;

     mass_part += _grad_test[_i][_qp](0) * _rho[0] * _u_vel_dot_du[_qp] * _phi[_j][_qp];

     Real mom_part = 0.0;

      mom_part += _grad_test[_i][_qp](0) * _rho[0] * ( vel * _grad_phi[_j][_qp] + _phi[_j][_qp] * _grad_u_vel[_qp](0) )
                 +_grad_test[_i][_qp](1) * _rho[1] * ( _phi[_j][_qp] * _grad_v_vel[_qp](0) )
		 +_grad_test[_i][_qp](2) * _rho[2] * ( _phi[_j][_qp] * _grad_w_vel[_qp](0) );
 
    Real diff_part = 0.0;
      
      diff_part -= _mu[0] * _grad_test[_i][_qp](0) * _second_phi[_j][_qp](_vert,_vert);

      diff_part -= _mu[0] * _coeff_param * _grad_test[_i][_qp](0) * _second_phi[_j][_qp](0,0);

      if(_vert > 1)
	  diff_part -= _mu[0] * _coeff_param * _grad_test[_i][_qp](0) * _second_phi[_j][_qp](1,1);


    return _tau_PSPG * ( mass_part + mom_part + diff_part);


    }else if (jvar == _v_vel_var_number){

	   RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);
       
	Real velmag2 = vel.norm_sq();
        
	Real _hele = 0.0;

	Real _muK = 0.0;

  if(_meshsize.size() == 1 ){
     _hele = _meshsize[0]; _muK = _mu[0];
  }else if(_meshsize.size() == 2){
     _hele = sqrt( _meshsize[0]*_meshsize[1] ); _muK = sqrt( _mu[0]*_mu[1] );
  }else{
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);
     _muK = std::pow(_mu[0]*_mu[1]*_mu[2],0.3333);
  }
  Real _tau_PSPG = _hele*_hele/_muK/4.0;

	Real mass_part = 0.0;

        mass_part += _grad_test[_i][_qp](1) * _rho[1] * _v_vel_dot_du[_qp] * _phi[_j][_qp];

  	Real mom_part = 0.0;
	Real diff_part = 0.0;
	
	mom_part += _grad_test[_i][_qp](0) * _rho[0] * ( _phi[_j][_qp] * _grad_u_vel[_qp](1) )
                      +_grad_test[_i][_qp](1) * _rho[1] * ( vel * _grad_phi[_j][_qp] + _phi[_j][_qp] * _grad_v_vel[_qp](1) )
		      +_grad_test[_i][_qp](2) * _rho[2] * ( _phi[_j][_qp] * _grad_w_vel[_qp](1) );
	   
	diff_part -= _mu[1] * _grad_test[_i][_qp](1) * _second_phi[_j][_qp](_vert,_vert);
	  
	 diff_part -= _mu[1] * _coeff_param * _grad_test[_i][_qp](1) * _second_phi[_j][_qp](0,0);
    
      if(_vert > 1)
          diff_part -= _mu[1] * _coeff_param * _grad_test[_i][_qp](1) * _second_phi[_j][_qp](1,1);




  	return _tau_PSPG * ( mass_part + mom_part + diff_part );
	       
	
    }else if(jvar == _w_vel_var_number){

	if(_vert > 1){

	RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);
       
 
	Real velmag2 = vel.norm_sq();

        Real _hele = 0.0;

	Real _muK = 0.0;

  if(_meshsize.size() == 1 ){
     _hele = _meshsize[0]; _muK = _mu[0];
  }else if(_meshsize.size() == 2){
     _hele = sqrt( _meshsize[0]*_meshsize[1] ); _muK = sqrt( _mu[0]*_mu[1] );
  }else{
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);
     _muK = std::pow(_mu[0]*_mu[1]*_mu[2],0.3333);
  }
  Real _tau_PSPG = _hele*_hele/_muK/4.0;


        Real mass_part = 0.0;
	Real mom_part  = 0.0;
	Real diff_part = 0.0;

	mass_part += _grad_test[_i][_qp](2) * _rho[2] * _w_vel_dot_du[_qp] * _phi[_j][_qp];
	
	mom_part += _grad_test[_i][_qp](0) * _rho[0] * ( _phi[_j][_qp] * _grad_u_vel[_qp](2) )
                   +_grad_test[_i][_qp](1) * _rho[1] * ( _phi[_j][_qp] * _grad_v_vel[_qp](2) )
                   +_grad_test[_i][_qp](2) * _rho[2] * ( vel * _grad_phi[_j][_qp] + _phi[_j][_qp] * _grad_w_vel[_qp](2) );		
	diff_part -= _mu[2] * _grad_test[_i][_qp](2) * _second_phi[_j][_qp](_vert,_vert);
	 diff_part -= _mu[2] * _coeff_param * _grad_test[_i][_qp](2) * _second_phi[_j][_qp](0,0);
    
      if(_vert > 1)
          diff_part -= _mu[2] * _coeff_param * _grad_test[_i][_qp](2) * _second_phi[_j][_qp](1,1);


	
	return _tau_PSPG * ( mass_part + mom_part + diff_part );
	}else
		return 0.0;
	 
    }else
	return 0.0;

}

