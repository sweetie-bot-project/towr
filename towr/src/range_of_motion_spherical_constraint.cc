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

#include <towr/constraints/range_of_motion_spherical_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

RangeOfMotionSphericalConstraint::RangeOfMotionSphericalConstraint (const AlignedBox3d& box,
                                                  const Sphere3d& sphere,
                                                  double T, double dt,
                                                  const EE& ee,
                                                  const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "rangeofmotionsherical-" + std::to_string(ee))
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_motion_    = spline_holder.ee_motion_.at(ee);

  ee_bounding_box_B_          = box;
  ee_bounding_sphere_B_       = sphere;
  ee_ = ee;

  SetRows(GetNumberOfNodes()*(k3D+1));
}

int
RangeOfMotionSphericalConstraint::GetRow (int node, int dim) const
{
  return node*(k3D+1) + dim;
}

void
RangeOfMotionSphericalConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector3d base_W  = base_linear_->GetPoint(t).p();
  Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  Vector3d vector_base_to_ee_B = b_R_w*(pos_ee_W - base_W);

  g.middleRows(GetRow(k, X), k3D) = vector_base_to_ee_B;
  g[GetRow(k, k3D)] = (vector_base_to_ee_B - ee_bounding_sphere_B_.center()).norm();
}

void
RangeOfMotionSphericalConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<k3D; ++dim) {
    bounds.at(GetRow(k,dim)) = ifopt::Bounds(ee_bounding_box_B_.min()[dim], ee_bounding_box_B_.max()[dim]);
  }
  bounds.at(GetRow(k,k3D)) = ifopt::Bounds(0.0, ee_bounding_sphere_B_.radius());
}

void
RangeOfMotionSphericalConstraint::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  int row_start = GetRow(k,X);

  Vector3d base_W  = base_linear_->GetPoint(t).p();
  Vector3d ee_W = ee_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  Vector3d vector_base_to_ee_B = b_R_w*(ee_W - base_W);
  Vector3d vector_center_to_ee_B = vector_base_to_ee_B - ee_bounding_sphere_B_.center();
  double distance_center_to_ee = vector_center_to_ee_B.norm();

  NodeSpline::Jacobian grad;
  if (var_set == id::base_lin_nodes) {
	grad = -1*b_R_w*base_linear_->GetJacobianWrtNodes(t, kPos);
  }
  else if (var_set == id::base_ang_nodes) {
    Vector3d vector_base_to_ee_W = ee_W - base_W;
	grad = base_angular_.DerivOfRotVecMult(t, vector_base_to_ee_W, true);
  }
  else if (var_set == id::EEMotionNodes(ee_)) {
    grad = b_R_w*ee_motion_->GetJacobianWrtNodes(t,kPos);
  }
  else if (var_set == id::EESchedule(ee_)) {
    grad = b_R_w*ee_motion_->GetJacobianOfPosWrtDurations(t);
  }
  else return;
 
  jac.middleRows(row_start, k3D) = grad;
  jac.middleRows(row_start+k3D, 1) = (1.0/distance_center_to_ee * vector_center_to_ee_B.transpose()).sparseView() * grad; 
}

} /* namespace xpp */
