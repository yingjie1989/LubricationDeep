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

#ifndef COUPLEDNEUMANNGRADIENTBC_H
#define COUPLEDNEUMANNGRADIENTBC_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannGradientBC;

template<>
InputParameters validParams<CoupledNeumannGradientBC>();

class CoupledNeumannGradientBC : public IntegratedBC
{
public:


 CoupledNeumannGradientBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 Real _Re;
 Real _tol;
 int _dim;
      const VariableValue & _some_var_val;
   const VariableGradient & _some_var_grad;

};

#endif //COUPLEDNEUMANNGRADIENTBC_H

