/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "InNSSUPGTimeDerivative.h"

template<>
InputParameters validParams<InNSSUPGTimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  //Coupled Variables

  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addCoupledVar("u_mesh", 0, "x-meshvelocity"); // only required in 3D
  params.addCoupledVar("v_mesh", 0, "y-meshvelocity"); // only required in 3D
  params.addCoupledVar("w_mesh", 0, "z-meshvelocity"); // only required in 3D


 // Required parameters
  params.addParam<Real>("timeStep", 0.01,"timeStep");
  params.addRequiredParam<Real>("mu", "dynamic viscosity");
  params.addParam<Real>("rho", "density");
  params.addParam<Real>("visco_coeff",1.0, "viscous coefficient");
  params.addRequiredParam<std::vector<Real> >("hele", "Coefficient of perturbation");
  params.addParam<bool>("SUPGold", false, "Use old velocity profile to stabilize");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addRequiredParam<unsigned>("vert", "vertical diffusion direction");
 
  return params;
}

InNSSUPGTimeDerivative::InNSSUPGTimeDerivative(const InputParameters & parameters) :
    TimeKernel(parameters),
   // Coupled variables
   _u_vel(coupledValue("u")),
   _v_vel(coupledValue("v")),
   _w_vel(coupledValue("w")),

   _u_vel_old(coupledValueOld("u")),
   _v_vel_old(coupledValueOld("v")),
   _w_vel_old(coupledValueOld("w")),

  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),

  _u_vel_mesh(coupledValue("u_mesh")),
  _v_vel_mesh(coupledValue("v_mesh")),
  _w_vel_mesh(coupledValue("w_mesh")),

  
  // Required parameters
  _DeltaT(getParam<Real>("timeStep")),
  _mu(getParam<Real>("mu")),
  _visco_coeff(getParam<Real>("visco_coeff")),
  _rho(getParam<Real>("rho")),
  _meshsize(getParam<std::vector<Real> >("hele")),
  _ifOld(getParam<bool>("SUPGold")),
  _component(getParam<unsigned>("component")),
  _vert(getParam<unsigned>("vert"))

{
}

Real
InNSSUPGTimeDerivative::computeQpResidual()
{
  // See "Component SUPG contributions" section of notes for details.

   //Real _u_vel_mesh = (_dispx_val[_qp] - _dispx_old[_qp])/_DeltaT;
   //Real _v_vel_mesh = (_dispy_val[_qp] - _dispy_old[_qp])/_DeltaT;
   //Real _w_vel_mesh = (_dispz_val[_qp] - _dispz_old[_qp])/_DeltaT;

   RealVectorValue vel;
   
   if(_ifOld){
	vel(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
	vel(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
	vel(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
   }else{
	vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
	vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
	vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];
   }
  // Velocity vector magnitude squared
  Real velmag2 = vel.norm_sq();

  // Velocity vector, dotted with the test function gradient

  Real _hele = 0.0;

  if(_meshsize.size() == 1 )
     _hele = _meshsize[0];
  else if(_meshsize.size() == 2)
     _hele = sqrt( _meshsize[0]*_meshsize[1] );
  else
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

  Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

  Real U_grad_phi = _rho * vel * _grad_test[_i][_qp];

  Real mass_part = _rho * _u_dot[_qp];

  Real result = _tau_SUPG * U_grad_phi * mass_part;

  return result;
}

Real
InNSSUPGTimeDerivative::computeQpJacobian()
{

   //Real _u_vel_mesh = (_dispx_val[_qp] - _dispx_old[_qp])/_DeltaT;
   //Real _v_vel_mesh = (_dispy_val[_qp] - _dispy_old[_qp])/_DeltaT;
   //Real _w_vel_mesh = (_dispz_val[_qp] - _dispz_old[_qp])/_DeltaT;

    RealVectorValue vel;
  
    if(_ifOld){
        vel(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
        vel(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
        vel(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
   }else{
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];
   }


  Real velmag2 = vel.norm_sq();
   
  Real _hele = 0.0;

  if(_meshsize.size() == 1 )
     _hele = _meshsize[0];
  else if(_meshsize.size() == 2)
     _hele = sqrt( _meshsize[0]*_meshsize[1] );
  else
     _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

   Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

   Real U_grad_phi = _rho * vel * _grad_test[_i][_qp];
 
   Real mass_part = _rho * _u_dot[_qp];

   Real mass_part_dev = _rho * _phi[_j][_qp]*_du_dot_du[_qp];
  
   Real part1 = _tau_SUPG * U_grad_phi * mass_part_dev;
   
   Real part2 = 0.0;
 
   if(!_ifOld){
	Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](_component) * _phi[_j][_qp];

	part2 += _tau_SUPG * U_grad_phi_dev * mass_part;
    }
   return part1 + part2;

}

Real
InNSSUPGTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{

    if(_ifOld)
	return 0.0;

    if (jvar == _u_vel_var_number){

          //Real _u_vel_mesh = (_dispx_val[_qp] - _dispx_old[_qp])/_DeltaT;
   	  //Real _v_vel_mesh = (_dispy_val[_qp] - _dispy_old[_qp])/_DeltaT;
          //Real _w_vel_mesh = (_dispz_val[_qp] - _dispz_old[_qp])/_DeltaT;

          RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);  	

	Real velmag2 = vel.norm_sq();

	 Real _hele = 0.0;

     if(_meshsize.size() == 1 )
     	_hele = _meshsize[0];
     else if(_meshsize.size() == 2)
      	_hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
     	_hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

   	Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);
	
   	Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](0) * _phi[_j][_qp];
        Real mass_part = _rho * _u_dot[_qp];

        Real part = _tau_SUPG * U_grad_phi_dev * mass_part;

        return part;

    }else if (jvar == _v_vel_var_number){

          //Real _u_vel_mesh = (_dispx_val[_qp] - _dispx_old[_qp])/_DeltaT;
   	  //Real _v_vel_mesh = (_dispy_val[_qp] - _dispy_old[_qp])/_DeltaT;
   	 //Real _w_vel_mesh = (_dispz_val[_qp] - _dispz_old[_qp])/_DeltaT;

	RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);

	Real velmag2 = vel.norm_sq();
        
	Real _hele = 0.0;
	if(_meshsize.size() == 1 )
           _hele = _meshsize[0];
     	else if(_meshsize.size() == 2)
           _hele = sqrt( _meshsize[0]*_meshsize[1] );
     	else
           _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

	Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);
        
        Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](1) * _phi[_j][_qp];
        Real mass_part = _rho * _u_dot[_qp];

        Real part = _tau_SUPG * U_grad_phi_dev * mass_part;

        return part;
	

    }else if (jvar == _w_vel_var_number){

          //Real _u_vel_mesh = (_dispx_val[_qp] - _dispx_old[_qp])/_DeltaT;
   	  //Real _v_vel_mesh = (_dispy_val[_qp] - _dispy_old[_qp])/_DeltaT;
   	  //Real _w_vel_mesh = (_dispz_val[_qp] - _dispz_old[_qp])/_DeltaT;

   RealVectorValue vel(_u_vel[_qp] - _u_vel_mesh[_qp], _v_vel[_qp] - _v_vel_mesh[_qp], _w_vel[_qp] - _w_vel_mesh[_qp]);

	Real velmag2 = vel.norm_sq();
        
	Real _hele = 0.0;
	if(_meshsize.size() == 1 )
          _hele = _meshsize[0];
        else if(_meshsize.size() == 2)
          _hele = sqrt( _meshsize[0]*_meshsize[1] );
        else
          _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

	Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

        Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](2) * _phi[_j][_qp];
        Real mass_part = _rho * _u_dot[_qp];

        Real part = _tau_SUPG * U_grad_phi_dev * mass_part;

        return part;
    }

    else 
	return 0.0;

}

