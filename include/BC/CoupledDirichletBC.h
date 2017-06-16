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

#ifndef COUPLEDDIRICHLETBC_H
#define COUPLEDDIRICHLETBC_H

#include "NodalBC.h"


class CoupledDirichletBC;

template<>
InputParameters validParams<CoupledDirichletBC>();
class CoupledDirichletBC : public NodalBC
{
public:CoupledDirichletBC(const InputParameters & parameters);

protected:

virtual Real computeQpResidual() override;

private:

Real _alpha;

const VariableValue & _some_var_val;
};


#endif //COUPLEDDIRICHLETBC_H
