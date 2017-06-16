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

#ifndef MODIFYREACTION_H
#define MODIFYREACTION_H

#include "Kernel.h"

// Forward Declaration
class ModifyReaction;

template<>
InputParameters validParams<ModifyReaction>();

class ModifyReaction : public Kernel
{
public:
  ModifyReaction(const InputParameters & parameters);

  Real _timeStep;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

};
#endif //MODIFYREACTION_H
