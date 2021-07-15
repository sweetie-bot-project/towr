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

#ifndef TOWR_NLP_FACTORY_PLANAR_H_
#define TOWR_NLP_FACTORY_PLANAR_H_

#include <towr/nlp_formulation_base.h>

namespace towr {

/**
 *
 * @brief Planar NLP formulation.
 *
 */
class NlpFormulationPlanar : public NlpFormulationBase {
public:
  NlpFormulationPlanar ();
  virtual ~NlpFormulationPlanar () = default;

  /**
   * @brief The ifopt variable sets that will be optimized over.
   * @param[in/out] builds fully-constructed splines from the variables.
   */
  VariablePtrVec GetVariableSets(SplineHolder& spline_holder);

  /**
   * @brief The ifopt constraints that enforce feasible motions.
   * @param[in] uses the fully-constructed splines for initialization of constraints.
   */
  ContraintPtrVec GetConstraints(const SplineHolder& spline_holder) const;

  /** @brief The ifopt costs to tune the motion. */
  ContraintPtrVec GetCosts(const SplineHolder& s) const;

private:
  // variables
  std::vector<NodesVariables::Ptr> MakeBaseVariables() const;
  std::vector<NodesVariablesPhaseBased::Ptr> MakeEndeffectorVariables() const;
  std::vector<NodesVariablesTimeBased::Ptr> MakeForceVariables() const;
  std::vector<PhaseDurations::Ptr> MakeContactScheduleVariables () const;

  // constraints
  ContraintPtrVec GetConstraint(Parameters::ConstraintName name,
                                const SplineHolder& splines) const;
  ContraintPtrVec MakeDynamicConstraint(const SplineHolder& s) const;
  ContraintPtrVec MakeRangeOfMotionConstraint(const SplineHolder& s) const;
  ContraintPtrVec MakeRangeOfMotionSphericalConstraint (const SplineHolder& s) const;
  ContraintPtrVec MakeTotalTimeConstraint() const;
  ContraintPtrVec MakeTerrainConstraint() const;
  ContraintPtrVec MakeForceConstraint() const;
  ContraintPtrVec MakeSwingConstraint() const;

  // costs
  CostPtrVec GetCost(const Parameters::CostName& id, double weight, const SplineHolder& s) const;
  CostPtrVec MakeForcesCost(double weight) const;
  CostPtrVec MakeEEMotionCost(double weight) const;
  CostPtrVec MakeBaseLinMotionCost(double weight) const;
  CostPtrVec MakeBaseAngMotionCost(double weight) const;
  CostPtrVec MakeEEAccCost(double weight, const SplineHolder& s) const;
};

} /* namespace towr */

#endif /* TOWR_NLP_FACTORY_H_ */
