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

#ifndef TOWR_MODELS_SIMPLE_KINEMATIC_MODEL_H_
#define TOWR_MODELS_SIMPLE_KINEMATIC_MODEL_H_

#include <towr/models/kinematic_model.h>

namespace towr {

/**
 * @brief Simplified kinematic model in form of bounding box.
 *
 * The dimentions of bounding boxes are the same for all end effectors.
 * Nominal pose is always in the senter of bounding box.
 *
 * @ingroup Robots
 */
class SimpleKinematicModel : public KinematicModel {
protected:
  // partial initaliztion for use in derived classes
  SimpleKinematicModel(int n_ee) {
    nominal_stance_.resize(n_ee); 
  }

public:
  // fully initialized Kinematic model
  SimpleKinematicModel(const EEPos& nominal_stance, const Vector3d& max_dev_from_nominal) :
	  nominal_stance_(nominal_stance), max_dev_from_nominal_(max_dev_from_nominal) 
  {}

  virtual ~SimpleKinematicModel () = default;

  // additional getter for backward compatibility
  Vector3d GetMaximumDeviationFromNominal() const {
    return max_dev_from_nominal_;
  }

  // kinematic model interface
  Vector3d GetNominalStanceInBase(EE ee) const { 
    return nominal_stance_.at(ee); 
  }

  AlignedBox3d GetBoundingBox(EE ee) const { 
    return AlignedBox3d(nominal_stance_.at(ee) - max_dev_from_nominal_, nominal_stance_.at(ee) + max_dev_from_nominal_); 
  }

  virtual Sphere3d GetBoundingBall(EE ee) const {
    return Sphere3d(Vector3d::Zero(), 0.0); 
  }

  virtual int GetNumberOfEndeffectors() const {
    return nominal_stance_.size();
  }

  protected:
    EEPos nominal_stance_;
    Vector3d max_dev_from_nominal_;
};

} /* namespace towr */

#endif /* TOWR_MODELS_SIMPLE_KINEMATIC_MODEL_H_ */
