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

#include <iostream>

#include <towr/constraints/range_of_motion_planar_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

RangeOfMotionPlanarConstraint::RangeOfMotionPlanarConstraint (const AlignedBox3d& box,
                                                  double T, double dt,
                                                  const EE& ee,
                                                  const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "rangeofmotion2d-" + std::to_string(ee))
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = spline_holder.base_angular_;
  ee_motion_    = spline_holder.ee_motion_.at(ee);

  ee_bounding_box_B_ = box;
  ee_ = ee;

  SetRows(GetNumberOfNodes()*k2D);
}

int
RangeOfMotionPlanarConstraint::GetRow (int node, int dim) const
{
  return node*k2D + dim;
}

void
RangeOfMotionPlanarConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector2d base_W  = base_linear_->GetPoint(t).p().head<k2D>();
  Vector2d pos_ee_W = ee_motion_->GetPoint(t).p().head<k2D>();
  double yaw = base_angular_->GetPoint(t).p().z();
  Eigen::Matrix2d b_R_w;
  b_R_w <<     cos(yaw),     sin(yaw),
		     - sin(yaw),     cos(yaw);

  Vector2d vector_base_to_ee_W = pos_ee_W - base_W;
  Vector2d vector_base_to_ee_B = b_R_w*(vector_base_to_ee_W);

  g.middleRows(GetRow(k, X), k2D) = vector_base_to_ee_B;
}

void
RangeOfMotionPlanarConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<k2D; ++dim) {
    ifopt::Bounds b;
    b.upper_ = ee_bounding_box_B_.max()[dim];
    b.lower_ = ee_bounding_box_B_.min()[dim];
    bounds.at(GetRow(k,dim)) = b;
  }
}

void
RangeOfMotionPlanarConstraint::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  double yaw = base_angular_->GetPoint(t).p().z();
  
  int row_start = GetRow(k,X);

  if (var_set == id::base_lin_nodes) {
    Eigen::SparseMatrix<double, Eigen::RowMajor> b_mRS_w(2,2);
    b_mRS_w.insert(0,0) = -cos(yaw);
	b_mRS_w.insert(0,1) = -sin(yaw);
    b_mRS_w.insert(1,0) = sin(yaw);
	b_mRS_w.insert(1,1) = -cos(yaw);

    jac.middleRows(row_start, k2D) = b_mRS_w * base_linear_->GetJacobianWrtNodes(t, kPos).topRows<k2D>();
  }

  if (var_set == id::base_ang_nodes) {
    Vector2d base_W   = base_linear_->GetPoint(t).p().head<k2D>();
    Vector2d ee_pos_W = ee_motion_->GetPoint(t).p().head<k2D>();
    Vector2d r_W = ee_pos_W - base_W;

    Eigen::SparseMatrix<double, Eigen::ColMajor> deriv(2, 1);
	deriv.insert(0,0) = - sin(yaw)*r_W.x() + cos(yaw)*r_W.y();
	deriv.insert(1,0) = - cos(yaw)*r_W.x() - sin(yaw)*r_W.y();

    jac.middleRows(row_start, k2D) = deriv * base_angular_->GetJacobianWrtNodes(t, kPos).row(Z);
  }

  if (var_set == id::EEMotionNodes(ee_)) {
    Eigen::SparseMatrix<double, Eigen::RowMajor> b_RS_w(2,2);
    b_RS_w.insert(0,0) = cos(yaw);
	b_RS_w.insert(0,1) = sin(yaw);
    b_RS_w.insert(1,0) = -sin(yaw);
	b_RS_w.insert(1,1) = cos(yaw);

    jac.middleRows(row_start, k2D) = b_RS_w*ee_motion_->GetJacobianWrtNodes(t,kPos).topRows<k2D>();
  }

  if (var_set == id::EESchedule(ee_)) {
    Eigen::SparseMatrix<double, Eigen::RowMajor> b_RS_w(2,2);
    b_RS_w.insert(0,0) = cos(yaw);
	b_RS_w.insert(0,1) = sin(yaw);
    b_RS_w.insert(1,0) = -sin(yaw);
	b_RS_w.insert(1,1) = cos(yaw);

    jac.middleRows(row_start, k2D) = b_RS_w*ee_motion_->GetJacobianOfPosWrtDurations(t).topRows<k2D>();
  }
}

} /* namespace xpp */

