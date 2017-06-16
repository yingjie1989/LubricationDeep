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

#ifndef COUPLEDNEUMANNGRADIENTBCMODIFY_H
#define COUPLEDNEUMANNGRADIENTBCMODIFY_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannGradientBCModify;

template<>
InputParameters validParams<CoupledNeumannGradientBCModify>();

class CoupledNeumannGradientBCModify : public IntegratedBC
{
public:


 CoupledNeumannGradientBCModify(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 	Real _Re;
 	Real _tol;
 	int _dim;
 	int _dirt;
 	
	const VariableValue & _some_var_val;
 	const VariableGradient & _c_grad;
        const VariableGradient & _h_grad;

};

#endif //COUPLEDNEUMANNGRADIENTBCMODIFY_H

