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


//include header of my application
#include "lubrication.h"

#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "MooseSyntax.h"
#include "NavierStokesApp.h"
//#include "MooseTestApp.h"


//modified BC

#include "CoupledNeumannGradientBC.h"
#include "CoupledNeumannGradientBCModify.h"
#include "CoupledNeumannCurveBC.h"
#include "CoupledNeumannCurveBCModify.h"
#include "CoupledNeumannBC.h"
#include "CoupledNeumannStressBC.h"
#include "CoupledDirichletBC.h"
#include "CoupledDirichletGradientBC.h"

//modified IC
#include "QuadraticIC.h"
#include "LinearIC.h"
#include "FingerIC.h"
#include "FingerCircleIC.h"
#include "tanhCircleIC.h"
#include "CrossIC.h"
#include "SmoothCircleIC.h"


//modified materials
#include "DerivativeParsedMaterial.h"
#include "SourceMaterial3D.h"
#include "SourceMaterial.h"
#include "CapilaryMaterial.h"
#include "CapilaryCircleMaterial.h"
#include "ComputeEnrichStrain.h"
#include "ComputeSurfaceStress.h"


//custom kernel
#include "Advection_c.h"
#include "Advection_v.h"
#include "Advection_h.h"
#include "Advection_hpp.h"
#include "Diffusion_D.h"
#include "DiffusionCoupled.h"
#include "ModifyReaction.h"
#include "ReactionCoupled.h"
#include "AdvectionSUPG.h"
#include "TimeDerivativeSUPG.h"
#include "TimeDerivative_D.h"
#include "ReactionSUPG.h"
#include "NSSelf.h"
#include "NSCoupled.h"
#include "NSDiffusion.h"
#include "InCompress.h"
#include "InCompressVelocity.h"
#include "InCompressStream.h"
#include "NSPressure.h"
#include "SourceDiff.h"
#include "SourceGradient.h"
#include "SourceManuFac.h"
#include "InNSMomentum.h"
#include "InNSSUPGMomentum.h"
#include "InNSSUPGTimeDerivative.h"
#include "InCompSUPG.h"
#include "InCompModify.h"
#include "INSMassModify.h"
#include "INSMomentumTractionFormNoDimen.h"


//auxkernel
#include "DeformAux.h"
#include "ResCheckAux.h"
#include "VelTransAux.h"
#include "SurfPresAux.h"
#include "VelTransGradAux.h"
#include "VelMeshCal.h"


//transfer
//#include "MultiAppDTKInterpolationTransfer.h"


template<>
InputParameters validParams<lubrication>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

lubrication::lubrication(InputParameters parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  lubrication::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  lubrication::associateSyntax(_syntax, _action_factory);
}

lubrication::~lubrication()
{
}

void
lubrication::registerObjects(Factory & factory)
{
//register objects from module phase_field
//kernels 

registerKernel(Diffusion_D);
registerKernel(Advection_c);
registerKernel(Advection_v);
registerKernel(Advection_h);
registerKernel(Advection_hpp);
registerKernel(AdvectionSUPG);
registerKernel(DiffusionCoupled);
registerKernel(ModifyReaction);
registerKernel(ReactionCoupled);
registerKernel(TimeDerivativeSUPG);
registerKernel(ReactionSUPG);
registerKernel(NSSelf);
registerKernel(NSCoupled);
registerKernel(NSDiffusion);
registerKernel(InCompress);
registerKernel(InCompressVelocity);
registerKernel(InCompressStream);
registerKernel(SourceDiff);
registerKernel(SourceGradient);
registerKernel(SourceManuFac);
registerKernel(NSPressure);
registerKernel(InNSMomentum);
registerKernel(InNSSUPGMomentum);
registerKernel(InNSSUPGTimeDerivative);
registerKernel(InCompSUPG);
registerKernel(InCompModify);
registerKernel(INSMassModify);
registerKernel(TimeDerivative_D);
//registerKernel(INSMomentumBaseALE);
registerKernel(INSMomentumTractionFormNoDimen);


//auxkernel
registerAux(DeformAux);
registerAux(ResCheckAux);
registerAux(VelTransAux);
registerAux(VelMeshCal);
registerAux(SurfPresAux);
registerAux(VelTransGradAux);


//IC
registerInitialCondition(QuadraticIC);
registerInitialCondition(LinearIC);
registerInitialCondition(FingerIC);
registerInitialCondition(FingerCircleIC);
registerInitialCondition(tanhCircleIC);
registerInitialCondition(SmoothCircleIC);
registerInitialCondition(CrossIC);

//BC
registerBoundaryCondition(CoupledNeumannBC);
registerBoundaryCondition(CoupledNeumannStressBC);
registerBoundaryCondition(CoupledNeumannGradientBC);
registerBoundaryCondition(CoupledNeumannGradientBCModify);
registerBoundaryCondition(CoupledNeumannCurveBC);
registerBoundaryCondition(CoupledNeumannCurveBCModify);
registerBoundaryCondition(CoupledDirichletBC);
registerBoundaryCondition(CoupledDirichletGradientBC);


//Materials
registerMaterial(SourceMaterial3D);
registerMaterial(SourceMaterial);
registerMaterial(CapilaryMaterial);
registerMaterial(CapilaryCircleMaterial);
registerMaterial(ComputeEnrichStrain);
registerMaterial(ComputeSurfaceStress);


//transfer
//registerObject(MultiAppDTKInterpolationTransfer);


NavierStokesApp::registerObjects(factory);
//MooseTestApp::registerObjects(factory);
}

void
lubrication::registerApps()
{
  registerApp(lubrication);
}

void
lubrication::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
    NavierStokesApp::associateSyntax(syntax, action_factory);
   // MooseTestApp::associateSyntax(syntax, action_factory);

}
