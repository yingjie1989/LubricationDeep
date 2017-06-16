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

#ifndef COUPLEDNEUMANNCURVEBC_H
#define COUPLEDNEUMANNCURVEBC_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannCurveBC;

template<>
InputParameters validParams<CoupledNeumannCurveBC>();

class CoupledNeumannCurveBC : public IntegratedBC
{
public:


 CoupledNeumannCurveBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 Real _Re;
 bool _ifOld;
 Real _tol;
 int _dim;
      const VariableValue & _some_var_val;
      const VariableValue & _p_coupled;
      const VariableValue & _p_coupled_old;
      const VariableGradient & _c_grad;
      const VariableGradient & _h_grad;
      const VariableSecond & _h_second;

};

#endif //COUPLEDNEUMANNGCURVEBC_H

