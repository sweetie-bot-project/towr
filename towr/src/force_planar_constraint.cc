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

#include <towr/constraints/force_planar_constraint.h>

#include <towr/variables/variable_names.h>

namespace towr {

ForcePlanarConstraint::ForcePlanarConstraint (double force_limit,
                                  EE ee)
    : ifopt::ConstraintSet(kSpecifyLater, "force2d-" + id::EEForceNodes(ee))
{
  fn_max_  = force_limit;
  ee_      = ee;
}

void
ForcePlanarConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  //TODO this function perform second stage of initizlization. In future this procedure should be unified for all constraints: 
  // either all initialization should happen in constructor, either all binding to variables should be happen hear.
  ee_force_  = x->GetComponent<NodesVariablesEEForceTimeBased>(id::EEForceNodes(ee_));
  // get list of optimized nodes
  pure_stance_force_node_ids_ = ee_force_->GetIndicesOfNonConstantNodes();
  // calculate number of constraints
  int constraint_count = pure_stance_force_node_ids_.size();
  SetRows(constraint_count);
}

Eigen::VectorXd
ForcePlanarConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row=0;
  const std::vector<towr::Node>& force_nodes = ee_force_->GetNodes();

  for (int f_node_id : pure_stance_force_node_ids_) {
	double f = force_nodes.at(f_node_id).p()[0];
    // unilateral force
    g(row++) = f; // >0 (unilateral forces)
  }

  return g;
}

ForcePlanarConstraint::VecBound
ForcePlanarConstraint::GetBounds () const
{
  VecBound bounds;

  for (int f_node_id : pure_stance_force_node_ids_) {
    bounds.push_back(ifopt::Bounds(0.0, fn_max_)); // unilateral forces
  }

  return bounds;
}

void
ForcePlanarConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{
  if (var_set == ee_force_->GetName()) {
    int row = 0;
    for (int f_node_id : pure_stance_force_node_ids_) {
      // unilateral force
	  // TODO use eye matrix

      int idx = ee_force_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, X));
	  assert(idx != NodesVariables::NodeValueNotOptimized);

	  jac.coeffRef(row++, idx) = 1.0;              // unilateral force
    }
  }
}


} /* namespace towr */
