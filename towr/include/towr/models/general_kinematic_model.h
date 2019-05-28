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

#ifndef TOWR_MODELS_GENERAL_KINEMATIC_MODEL_H_
#define TOWR_MODELS_GENERAL_KINEMATIC_MODEL_H_

#include <bitset>
#include <towr/models/kinematic_model.h>

namespace towr {

/**
 * @brief Complex kineamtic model with individual constraints for each leg.
 *
 * Movement of legs is restricted by bounding box and ball with given center.
 *
 * @ingroup Robots
 */

class GeneralKinematicModel : public KinematicModel {
public:
  static const size_t max_number_of_legs_ = 32;

public:
  GeneralKinematicModel(int n_ee);
  virtual ~GeneralKinematicModel () = default;

  /**
   * @brief Set end effector kinematic restrictions.
   * @param ee End effector index.
   * @param nominal_stance Nominal pose of end effector in base coordinates.
   * @param bounding_box Range restrictions in form of parallelepiped (Eigen::AlignedBox) in base frame.
   * @param bounding_sphere Range restrictions in form of sphere (Eigen::AlignedBox) in base frame.
   **/
  void configureEndEffector(EE ee, const Vector3d& nominal_stance, const AlignedBox3d& bounding_box, const Sphere3d& bounding_sphere); 

  /**
   * @brief Return true if all end effector were configured.
   * @return true if all end effectors were initializated.
   **/
  bool isConfigured() {
    int n_ee = nominal_stance_.size();
    return configured_legs_.count() == n_ee;
  }

  virtual Vector3d GetNominalStanceInBase(EE ee) const {
    return nominal_stance_.at(ee);
  }

  virtual AlignedBox3d GetBoundingBox(EE ee) const {
    return bounding_box_.at(ee);
  }

  virtual Sphere3d GetBoundingBall(EE ee) const {
    return bounding_sphere_.at(ee);
  }

  virtual int GetNumberOfEndeffectors() const {
    return nominal_stance_.size();
  }

protected:
  std::bitset<max_number_of_legs_> configured_legs_;
  std::vector<Vector3d> nominal_stance_;
  std::vector<AlignedBox3d> bounding_box_;
  std::vector<Sphere3d> bounding_sphere_;
};

} /* namespace towr */

#endif /* TOWR_MODELS_GENERAL_KINEMATIC_MODEL_H_ */
