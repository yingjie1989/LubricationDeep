/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "InNSSUPGMomentum.h"

template<>
InputParameters validParams<InNSSUPGMomentum>()
{
  InputParameters params = validParams<Kernel>();
  //Coupled Variables

  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addCoupledVar("p", 0, "pressure"); // pressure field
 
  params.addCoupledVar("u_mesh", 0, "x-meshvelocity"); // only required in 3D
  params.addCoupledVar("v_mesh", 0, "y-meshvelocity"); // only required in 3D
  params.addCoupledVar("w_mesh", 0, "z-meshvelocity"); // only required in 3D

 // Required parameters
  params.addParam<Real>("timeStep",1.0, "timeStep");
  params.addRequiredParam<Real>("mu", "dynamic viscosity");
  params.addParam<Real>("rho",1.0, "density");
  params.addParam<Real>("coeff_pressure",1.0, "pressure coefficient");
  params.addParam<Real>("coeff_param",0.0, "scale coefficient");
  params.addParam<Real>("grav", 0.0, "gravity");
  params.addParam<Real>("visco_coeff", 1.0, "viscous coefficient");
  params.addRequiredParam<std::vector<Real> >("hele", "Coefficient of perturbation");  
  params.addParam<bool>("SUPGold", false, "Use old velocity profile to stabilize");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addRequiredParam<unsigned>("vert", "vertical diffusion direction");
 
  return params;
}

InNSSUPGMomentum::InNSSUPGMomentum(const InputParameters & parameters) :
    Kernel(parameters),
   // Coupled variables
   _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second()),

   _u_vel(coupledValue("u")),
   _v_vel(coupledValue("v")),
   _w_vel(coupledValue("w")),
   _p(coupledValue("p")),

  _u_vel_mesh(coupledValue("u_mesh")),
  _v_vel_mesh(coupledValue("v_mesh")),
  _w_vel_mesh(coupledValue("w_mesh")),
   
   _u_vel_old(coupledValueOld("u")),
   _v_vel_old(coupledValueOld("v")),
   _w_vel_old(coupledValueOld("w")),
  
  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),
  _grad_p(coupledGradient("p")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _p_var_number(coupled("p")),


  // Required parameters
  _DeltaT(getParam<Real>("timeStep")),
  _mu(getParam<Real>("mu")),
  _visco_coeff(getParam<Real>("visco_coeff")),
  _grav(getParam<Real>("grav")),
  _rho(getParam<Real>("rho")),
  _meshsize(getParam<std::vector<Real> >("hele")),
  _ifOld(getParam<bool>("SUPGold")),
  _component(getParam<unsigned>("component")),
  _coeff(getParam<Real>("coeff_pressure")),
  _coeff_param(getParam<Real>("coeff_param")),
  _vert(getParam<unsigned>("vert"))
  
{
}

Real
InNSSUPGMomentum::computeQpResidual()
{
  // See "Component SUPG contributions" section of notes for details.

    RealVectorValue vel;
    RealVectorValue velOld;

   if(_ifOld){
        velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
        velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
        velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
   }
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];

  Real velmag2 = vel.norm_sq();
  Real U_grad_phi;
  
  if(_ifOld)
      U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
  else
      U_grad_phi = _rho * vel * _grad_test[_i][_qp];


  Real _hele = 0.0;
  if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);


  Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);
 
  Real convective_part = _rho * vel * _grad_u[_qp];

  Real pressure_part = _coeff * _grad_p[_qp](_component);

  Real diffusive_part = -_mu * _second_u[_qp](_vert,_vert);
 
  diffusive_part -= _mu * _coeff_param * _second_u[_qp](0,0);

  if(_vert > 1)
	diffusive_part -= _mu * _coeff_param * _second_u[_qp](1,1);


  Real result = _tau_SUPG * U_grad_phi * ( convective_part + diffusive_part + pressure_part + _grav );

  return result;
}

Real
InNSSUPGMomentum::computeQpJacobian()
{

    RealVectorValue vel;
    RealVectorValue velOld;

   if(_ifOld){
        velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
        velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
        velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
   }
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];
 
   Real velmag2 = vel.norm_sq();
 
   Real _hele = 0.0;
   if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

   Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

   Real U_grad_phi;
    
   if(_ifOld)
      U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
    else
      U_grad_phi = _rho * vel * _grad_test[_i][_qp];

   
   Real convective_part = _rho * ( (_u_vel[_qp] - _u_vel_mesh[_qp])*_grad_u[_qp](0) + (_v_vel[_qp] - _v_vel_mesh[_qp])*_grad_u[_qp](1) + (_w_vel[_qp] - _w_vel_mesh[_qp]) *_grad_u[_qp](2));
 
   //Real convective_part = _rho * vel * _grad_u[_qp];
   Real pressure_part = _coeff * _grad_p[_qp](_component);
   Real diffusive_part = -_mu * _second_u[_qp](_vert,_vert);

   diffusive_part -= _mu * _coeff_param * _second_u[_qp](0,0);

  if(_vert > 1)
        diffusive_part -= _mu * _coeff_param * _second_u[_qp](1,1);   

   Real convective_part_dev = _rho * ( (vel*_grad_phi[_j][_qp]) + _phi[_j][_qp]*_grad_u[_qp](_component) );
   Real diffusive_part_dev = -_mu * _second_phi[_j][_qp](_vert,_vert);  

   diffusive_part_dev -= _mu * _coeff_param * _second_phi[_j][_qp](0,0);

  if(_vert > 1)
        diffusive_part_dev -= _mu * _coeff_param * _second_phi[_j][_qp](1,1);

   Real part1 = _tau_SUPG * U_grad_phi * ( convective_part_dev + diffusive_part_dev );
   Real part2 = 0.0;

   if(!_ifOld){ 
     Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](_component) * _phi[_j][_qp];
     part2 += _tau_SUPG * U_grad_phi_dev * ( convective_part + diffusive_part + pressure_part + _grav ) ;   
   }

   return part1 + part2;

}

Real
InNSSUPGMomentum::computeQpOffDiagJacobian(unsigned int jvar)
{

    if (jvar == _u_vel_var_number){

	  RealVectorValue vel;
    	RealVectorValue velOld;

   	if(_ifOld){
        	velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
        	velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
        	velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
   	}
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];
        
	Real velmag2 = vel.norm_sq();

	Real _hele = 0.0;
	if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

   	Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);
	
	  Real U_grad_phi;
   
  	 if(_ifOld)
      		U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
    	else
      		U_grad_phi = _rho * vel * _grad_test[_i][_qp];
   	
   	//Real convective_part = _rho * (_u_vel[_qp]*_grad_u[_qp](0) + _v_vel[_qp]*_grad_u[_qp](1) + _w_vel[_qp]*_grad_u[_qp](2));	
	
	Real convective_part = _rho * vel * _grad_u[_qp];
	Real pressure_part = _coeff * _grad_p[_qp](_component);
        Real diffusive_part = -_mu * _second_u[_qp](_vert,_vert);

	 diffusive_part -= _mu * _coeff_param * _second_u[_qp](0,0);

  	if(_vert > 1)
        	diffusive_part -= _mu * _coeff_param * _second_u[_qp](1,1);


   	Real convective_part_dev = _rho * _phi[_j][_qp] * _grad_u[_qp](0);
		

	Real part1 = _tau_SUPG * U_grad_phi * convective_part_dev;

	Real part2 = 0.0;

	if(!_ifOld){
		Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](0) * _phi[_j][_qp];
		part2 += _tau_SUPG * U_grad_phi_dev * ( convective_part + diffusive_part + pressure_part + _grav);
	}

   	return part1 + part2;

    }else if (jvar == _v_vel_var_number){

	   RealVectorValue vel;
        RealVectorValue velOld;

        if(_ifOld){
                velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
                velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
                velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
        }
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];

	Real velmag2 = vel.norm_sq();

	Real _hele =0.0;
	if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

        Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);
        
	Real U_grad_phi;

         if(_ifOld)
                U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
        else
                U_grad_phi = _rho * vel * _grad_test[_i][_qp];
        
        //Real convective_part = _rho * (_u_vel[_qp]*_grad_u[_qp](0) + _v_vel[_qp]*_grad_u[_qp](1) + _w_vel[_qp]*_grad_u[_qp](2));
        Real convective_part = _rho * vel * _grad_u[_qp];
	Real pressure_part = _coeff * _grad_p[_qp](_component);
        Real diffusive_part = -_mu * _second_u[_qp](_vert,_vert);
	 diffusive_part -= _mu * _coeff_param * _second_u[_qp](0,0);

  	if(_vert > 1)
        	diffusive_part -= _mu * _coeff_param * _second_u[_qp](1,1);


        Real convective_part_dev = _rho * _phi[_j][_qp] * _grad_u[_qp](1);

        Real part1 = _tau_SUPG * U_grad_phi * convective_part_dev;

	 Real part2 = 0.0;

        if(!_ifOld){
        	Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](1) * _phi[_j][_qp];
		part2 += _tau_SUPG * U_grad_phi_dev * ( convective_part + diffusive_part + pressure_part + _grav );
	}


	return part1 + part2;
	

    }else if (jvar == _w_vel_var_number){

	   RealVectorValue vel;
        RealVectorValue velOld;

        if(_ifOld){
                velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
                velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
                velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
        }
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];

	Real velmag2 = vel.norm_sq();

	Real _hele = 0.0;
	if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

        Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

        Real U_grad_phi;

         if(_ifOld)
                U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
        else
                U_grad_phi = _rho * vel * _grad_test[_i][_qp];

        //Real convective_part = _rho * (_u_vel[_qp]*_grad_u[_qp](0) + _v_vel[_qp]*_grad_u[_qp](1) + _w_vel[_qp]*_grad_u[_qp](2));
	Real convective_part = _rho * vel * _grad_u[_qp];
	Real pressure_part = _coeff * _grad_p[_qp](_component);
        Real diffusive_part = -_mu * _second_u[_qp](_vert,_vert);        
	 diffusive_part -= _mu * _coeff_param * _second_u[_qp](0,0);

  	if(_vert > 1)
        	diffusive_part -= _mu * _coeff_param * _second_u[_qp](1,1);


	Real convective_part_dev = _rho * _phi[_j][_qp] * _grad_u[_qp](2);

        Real part1 = _tau_SUPG * U_grad_phi * convective_part_dev;
	Real part2 = 0.0;

        if(!_ifOld){
                Real U_grad_phi_dev = _rho * _grad_test[_i][_qp](2) * _phi[_j][_qp];
		part2 += _tau_SUPG * U_grad_phi_dev * ( convective_part + diffusive_part + pressure_part + _grav );
        }

	return part1 + part2;
    }
    else if(jvar == _p_var_number){

           RealVectorValue vel;
        RealVectorValue velOld;

        if(_ifOld){
                velOld(0) =  _u_vel_old[_qp] - _u_vel_mesh[_qp];
                velOld(1) =  _v_vel_old[_qp] - _v_vel_mesh[_qp];
                velOld(2) =  _w_vel_old[_qp] - _w_vel_mesh[_qp];
        }
        vel(0) = _u_vel[_qp] - _u_vel_mesh[_qp];
        vel(1) = _v_vel[_qp] - _v_vel_mesh[_qp];
        vel(2) = _w_vel[_qp] - _w_vel_mesh[_qp];

	Real velmag2 = vel.norm_sq();

	Real U_grad_phi;

         if(_ifOld)
                U_grad_phi = _rho * velOld * _grad_test[_i][_qp];
        else
                U_grad_phi = _rho * vel * _grad_test[_i][_qp];

        Real _hele = 0.0;
	if(_meshsize.size() == 1 )
        _hele = _meshsize[0];
     else if(_meshsize.size() == 2)
        _hele = sqrt( _meshsize[0]*_meshsize[1] );
     else
        _hele = std::pow(_meshsize[0]*_meshsize[1]*_meshsize[2],0.3333);

  	Real _tau_SUPG = _visco_coeff * _hele/sqrt( 16*_mu*_mu/_hele/_hele + 4*_rho*_rho*velmag2);

	Real pressure_part_dev = _coeff * _grad_phi[_j][_qp](_component);
		
	Real part = _tau_SUPG * U_grad_phi * pressure_part_dev;

	return part;

    }else
	return 0.0;

}

