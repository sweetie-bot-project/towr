/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/constraints/dynamic_planar_constraint.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

DynamicZMPPlanarConstraint::DynamicZMPPlanarConstraint (const DynamicModel::Ptr& m,
                                      double T, double dt,
                                      const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "dynamic-zmp-2d")
{
  base_mass_ = m->m();
  g_ = - m->g();

  // link with up-to-date spline variables
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = spline_holder.base_angular_;
  ee_motion_    = spline_holder.ee_motion_;
  ee_force_    = spline_holder.ee_force_;

  // iterate over time and count leg in contact

  SetRows(GetNumberOfNodes()*k3D);
}

int
DynamicZMPPlanarConstraint::GetRow (int k, int dim) const
{
  return k3D*k + dim;
}

void
DynamicZMPPlanarConstraint::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  Eigen::Vector3d v;
  State s = base_linear_->GetPoint(t);

  // ZMP x,y coordinates multiplied on Mg for base
  v.head<k2D>() = base_mass_ * g_ * s.p().head<k2D>() - base_mass_ * s.p().z() * s.a().head<k2D>();
  // ZMP ma coefficient
  v.z() = base_mass_ * g_;

  for(int ee=0; ee<ee_motion_.size(); ++ee)
  {
    s = ee_motion_[ee]->GetPoint(t);

    // take in account end effector masses
    if (ee_masses_.size() != 0) {
	  v.head<k2D>() += ee_masses_[ee] * (g_ - s.a().z()) * s.p().head<k2D>();
	  v.head<k2D>() -= ee_masses_[ee] * s.p().z() * s.a().head<k2D>();
	  v.z() += ee_masses_[ee] * (g_ - s.a().z());
	}
    // reaction forces
	double f = ee_force_[ee]->GetPoint(t).p().x();
	v.head<k2D>() += f * s.p().head<k2D>();
	v.z() += f;
    
  }
  g.segment(GetRow(k,AX), k3D) = v;
}

void
DynamicZMPPlanarConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const
{
  for (int dim = 0; dim < k3D; ++dim)
    bounds.at(GetRow(k,dim)) = ifopt::BoundZero;
}

void
DynamicZMPPlanarConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set,
                                            Jacobian& jac) const
{
  int n = jac.cols();
  Jacobian jac_model(k3D,n);


  // sensitivity of dynamic constraint w.r.t base variables.
  if (var_set == id::base_lin_nodes) {
    State s = base_linear_->GetPoint(t);

    Jacobian jac_base_lin_pos = base_linear_->GetJacobianWrtNodes(t,kPos);
    Jacobian jac_base_lin_acc = base_linear_->GetJacobianWrtNodes(t,kAcc);

    jac_model.topRows<k2D>() = base_mass_ * g_ * jac_base_lin_pos.topRows<k2D>() - base_mass_ * s.p().z() * jac_base_lin_acc.topRows<k2D>();
  }

  // sensitivity of dynamic constraint w.r.t. endeffector variables
  for (int ee=0; ee<ee_motion_.size(); ++ee) {
    if (var_set == id::EEForceNodes(ee)) {
      State s = ee_motion_[ee]->GetPoint(t);

      Jacobian jac_ee_force = ee_force_.at(ee)->GetJacobianWrtNodes(t,kPos);
      jac_model.row(X) = s.p().x() * jac_ee_force;
      jac_model.row(Y) = s.p().y() * jac_ee_force;
      jac_model.row(Z) = jac_ee_force;
	}

    if (var_set == id::EEMotionNodes(ee)) {
	  State s = ee_motion_.at(ee)->GetPoint(t);
	  double f = ee_force_.at(ee)->GetPoint(t).p().x();

      Jacobian jac_ee_pos = ee_motion_.at(ee)->GetJacobianWrtNodes(t,kPos);
      jac_model.topRows<k2D>() = f * jac_ee_pos.topRows<k2D>();

	  if (ee_masses_.size() != 0) {
        Jacobian jac_ee_acc = ee_motion_.at(ee)->GetJacobianWrtNodes(t,kAcc);

	    jac_model.row(X) -= (ee_masses_[ee] * s.p().x()) * jac_ee_acc.row(Z);
	    jac_model.row(Y) -= (ee_masses_[ee] * s.p().y()) * jac_ee_acc.row(Z);
		jac_model.row(X) -= (ee_masses_[ee] * s.a().x()) * jac_ee_pos.row(Z);
		jac_model.row(Y) -= (ee_masses_[ee] * s.a().y()) * jac_ee_pos.row(Z);
        jac_model.topRows<k2D>() += (ee_masses_[ee] * (g_ - s.a().z())) * jac_ee_pos.topRows<k2D>();
	    jac_model.topRows<k2D>() -= (ee_masses_[ee] * s.p().z()) * jac_ee_acc.topRows<k2D>();
	    jac_model.row(Z) = (-ee_masses_[ee]) * jac_ee_acc.row(Z);
	  }
    }
  }

  jac.middleRows(GetRow(k,AX), k3D) = jac_model;
}

} /* namespace towr */
