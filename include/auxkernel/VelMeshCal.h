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

#ifndef VELMESHCAL_H
#define VELMESHCAL_H

#include "AuxKernel.h"

class VelMeshCal;

template<>
InputParameters validParams<VelMeshCal>();

class VelMeshCal : public AuxKernel
{
public:

 VelMeshCal(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled;
  const VariableValue & _coupledOld;

  Real _DeltaT;

};


#endif //VELMESHCal_H
