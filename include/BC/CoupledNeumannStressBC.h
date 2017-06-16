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

#ifndef COUPLEDNEUMANNSTRESSBC_H
#define COUPLEDNEUMANNSTRESSBC_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannStressBC;

template<>
InputParameters validParams<CoupledNeumannStressBC>();

class CoupledNeumannStressBC : public IntegratedBC
{
public:


 CoupledNeumannStressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:

      const MaterialProperty<RealVectorValue> & _surface_stress;
      unsigned _dim;
      unsigned _component;
      const VariableValue & _some_var_val;
      const VariableGradient & _c_grad;
      const VariableGradient & _grad_h;

};

#endif //COUPLEDNEUMANNSTRESSBC_H

