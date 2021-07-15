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

#ifndef TOWR_CONSTRAINTS_RANGE_OF_MOTION_SPHERICAL_PLANAR_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_RANGE_OF_MOTION_SPHERICAL_PLANAR_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/sphere.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include "time_discretization_constraint.h"

namespace towr {

/** @brief Constrains an endeffector to lie in intersection of box and sphere.
  *
  * These constraints are necessary to avoid configurations
  * that are outside the kinematic reach of the robot. The constraint
  * is defined by Cartesian estimates of the reachability of each endeffector.
  *
  * This constraint calculates the position of of the contact expressed in the
  * current CoM frame and constrains it to lie innside intersection of box and sphere.
  *
  * @ingroup Constraints
  */
class RangeOfMotionSphericalPlanarConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector2d = Eigen::Vector2d;
  using Vector3d = Eigen::Vector3d;
  using AlignedBox3d = Eigen::AlignedBox3d;

  /**
   * @brief Constructs a constraint instance.
   * @param box   The bounding box for end effector.
   * @param sphere The bounding sphere.
   * @param T   The total duration of the optimization.
   * @param dt  the discretization intervall at which to enforce constraints.
   * @param ee            The endeffector for which to constrain the range.
   * @param spline_holder Pointer to the current variables.
   */
  RangeOfMotionSphericalPlanarConstraint(const AlignedBox3d& box,
                          const Sphere3d& sphere,
                          double T, double dt,
                          const EE& ee,
                          const SplineHolder& spline_holder);
  virtual ~RangeOfMotionSphericalPlanarConstraint() = default;

private:
  NodeSpline::Ptr base_linear_;     ///< the linear position of the base.
  NodeSpline::Ptr base_angular_;    ///< the orientation of the base.
  NodeSpline::Ptr ee_motion_;       ///< the linear position of the endeffectors.

  AlignedBox3d ee_bounding_box_B_; ///< bounding box for end effector position in base frame
  Sphere3d ee_bounding_sphere_B_; ///< bounding sphere for end effector position in base frame
  EE ee_; ///< index of end effector assotiated with constraint

  // see TimeDiscretizationPlanarConstraint for documentation
  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_RANGE_OF_MOTION_SPHERICAL_PLANAR_CONSTRAINT_H_ */