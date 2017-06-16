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

#include "TimeDerivative_D.h"
#include "Assembly.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<TimeDerivative_D>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addParam<Real>("coeff", 1.0, "viscous regularity coefficient");
  return params;
}

TimeDerivative_D::TimeDerivative_D(const InputParameters & parameters) :
    TimeKernel(parameters),
    _coeff(getParam<Real>("coeff")),
    _lumping(getParam<bool>("lumping"))
{
}

Real
TimeDerivative_D::computeQpResidual()
{
  return _coeff*_test[_i][_qp]*_u_dot[_qp];
}

Real
TimeDerivative_D::computeQpJacobian()
{
  return _coeff*_test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}

void
TimeDerivative_D::computeJacobian()
{
  if (_lumping)
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
  }
  else
    TimeKernel::computeJacobian();
}

