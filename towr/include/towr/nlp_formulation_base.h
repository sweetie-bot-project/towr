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

#ifndef TOWR_NLP_FACTORY_BASE_H_
#define TOWR_NLP_FACTORY_BASE_H_

#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

#include <towr/variables/spline_holder.h>
#include <towr/models/robot_model.h>
#include <towr/terrain/height_map.h>
#include <towr/parameters.h>

namespace towr {

/**
 * @defgroup Constraints
 * @brief Constraints of the trajectory optimization problem.
 *
 * These are the constraint sets that characterize legged locomotion.
 *
 * Folder: @ref include/towr/constraints
 */

/**
 * @defgroup Costs
 * @brief Costs of the trajectory optimization problem.
 *
 * These the the cost terms that prioritize certain solutions to the
 * legged locomotion problem.
 *
 * Folder: @ref include/towr/costs
 */

/**
 *
 * @brief A sample combination of variables, cost and constraints.
 *
 * This is _one_ example of how to combine the variables, constraints and costs
 * provided by this library. Additional variables or constraints can be added
 * to the NLP, or existing elements replaced to find a more powerful/general
 * formulation. This formulation was used to generate the motions described
 * in this paper: https://ieeexplore.ieee.org/document/8283570/
 */
class NlpFormulationBase {
public:
  using VariablePtrVec   = std::vector<ifopt::VariableSet::Ptr>;
  using ContraintPtr     = ifopt::ConstraintSet::Ptr;
  using ContraintPtrVec  = std::vector<ifopt::ConstraintSet::Ptr>;
  using CostPtrVec       = std::vector<ifopt::CostTerm::Ptr>;
  using EEPos            = std::vector<Eigen::Vector3d>;
  using Vector3d         = Eigen::Vector3d;

  virtual ~NlpFormulationBase () {};

  /**
   * @brief The ifopt variable sets that will be optimized over.
   * @param[in/out] builds fully-constructed splines from the variables.
   */
  virtual VariablePtrVec GetVariableSets(SplineHolder& spline_holder) = 0;

  /**
   * @brief The ifopt constraints that enforce feasible motions.
   * @param[in] uses the fully-constructed splines for initialization of constraints.
   */
  virtual ContraintPtrVec GetConstraints(const SplineHolder& spline_holder) const = 0;

  /** @brief The ifopt costs to tune the motion. */
  virtual ContraintPtrVec GetCosts(const SplineHolder& s) const = 0;

  BaseState initial_base_;
  BaseState final_base_;
  EEPos  initial_ee_W_;
  EEPos  final_ee_W_;
  RobotModel model_;
  HeightMap::Ptr terrain_;
  Parameters params_;
};

} /* namespace towr */

#endif /* TOWR_NLP_FACTORY_BASE_H_ */
