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

#ifndef TOWR_MODELS_KINEMATIC_MODEL_H_
#define TOWR_MODELS_KINEMATIC_MODEL_H_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <towr/variables/sphere.h>

namespace towr {

/**
 * @brief Contains all the robot specific kinematic parameters.
 *
 * This class is mainly used to formulate the @ref RangeOfMotionConstraint,
 * restricting each endeffector to stay inside it's kinematic range. 
 * The inposed restrictions include parallelepiped (bounding box) and ball 
 * with center in robot shoulder.
 *
 * @ingroup Robots
 */
class KinematicModel {
public:
  using Ptr      = std::shared_ptr<KinematicModel>;
  using EEPos    = std::vector<Eigen::Vector3d>;
  using Vector3d = Eigen::Vector3d;
  using AlignedBox3d    = Eigen::AlignedBox3d;
  using EE       = unsigned int;

  virtual ~KinematicModel () = default;

  /**
   * @brief  The xyz-position [m] of foot in default stance.
   * @param ee End effector index.
   * @returns The foot pose expressed in the base frame.
   */
  virtual Vector3d GetNominalStanceInBase(EE ee) const = 0;

  /**
   * @brief Get foot bounding box. 
   * @param ee End effector index.
   * @returns Pair of points which represent lower (first element) and upper (second element) bounds for end effector pose.
   */
  virtual AlignedBox3d GetBoundingBox(EE ee) const = 0;

  /**
   * @brief Get foot bounding sphere for the foot.
   * @param ee End effector index.
   * @return The center and raius of ball which restricts foot movements.
   */
  virtual Sphere3d GetBoundingSphere(EE ee) const = 0;

  /**
   * @returns returns the number of endeffectors of this robot.
   */
  virtual int GetNumberOfEndeffectors() const = 0;
};

} /* namespace towr */

#endif /* TOWR_MODELS_KINEMATIC_MODEL_H_ */
