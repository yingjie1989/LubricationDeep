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

#ifndef COUPLEDDIRICHLETGRADIENTBC_H
#define COUPLEDDIRICHLETGRADIENTBC_H

#include "NodalBC.h"


class CoupledDirichletGradientBC;

template<>
InputParameters validParams<CoupledDirichletGradientBC>();
class CoupledDirichletGradientBC : public NodalBC
{
public:CoupledDirichletGradientBC(const InputParameters & parameters);

protected:

virtual Real computeQpResidual() override;

private:

Real _alpha;
int  _component;
const VariableGradient & _some_var_grad;
};


#endif //COUPLEDDIRICHLETGRADIENTBC_H
