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

#ifndef COUPLEDNEUMANNCURVEBCMODIFY_H
#define COUPLEDNEUMANNCURVEBCMODIFY_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannCurveBCModify;

template<>
InputParameters validParams<CoupledNeumannCurveBCModify>();

class CoupledNeumannCurveBCModify : public IntegratedBC
{
public:


 CoupledNeumannCurveBCModify(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 Real _Re;
 bool _ifOld;
 Real _tol;
 unsigned _dim;
 unsigned _component;
      const VariableValue & _some_var_val;
      const VariableValue & _p_coupled;
      const VariableValue & _p_coupled_old;
      const VariableGradient & _c_grad;
      const VariableGradient & _grad_h;
      const VariableSecond & _h_second;

};

#endif //COUPLEDNEUMANNGCURVEBCMODIFY_H

