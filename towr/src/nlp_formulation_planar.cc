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

#include <towr/nlp_formulation_planar.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/phase_durations.h>

#include <towr/constraints/dynamic_planar_constraint.h>
#include <towr/constraints/force_planar_constraint.h>
#include <towr/constraints/range_of_motion_planar_constraint.h>
#include <towr/constraints/range_of_motion_spherical_planar_constraint.h>
#include <towr/constraints/swing_constraint.h>

#include <towr/costs/node_cost.h>
#include <towr/costs/accel_square_norm_cost.h>
#include <towr/variables/nodes_variables_partial.h>

#include <iostream>

namespace towr {

NlpFormulationPlanar::NlpFormulationPlanar ()
{
  using namespace std;
  cout << "\n";
  cout << "************************************************************\n";
  cout << " TOWR - Trajectory Optimization for Walking Robots (v1.4)\n";
  cout << "                \u00a9 Alexander W. Winkler\n";
  cout << "           https://github.com/ethz-adrl/towr\n";
  cout << "************************************************************";
  cout << "\n\n";
}

NlpFormulationPlanar::VariablePtrVec
NlpFormulationPlanar::GetVariableSets (SplineHolder& spline_holder)
{
  VariablePtrVec vars;

  auto base_motion = MakeBaseVariables();
  vars.insert(vars.end(), base_motion.begin(), base_motion.end());

  auto ee_motion = MakeEndeffectorVariables();
  vars.insert(vars.end(), ee_motion.begin(), ee_motion.end());

  auto ee_force = MakeForceVariables();
  vars.insert(vars.end(), ee_force.begin(), ee_force.end());

  auto contact_schedule = MakeContactScheduleVariables();

  // stores these readily constructed spline
  spline_holder = SplineHolder(base_motion.at(0), // linear
                               base_motion.at(1), // angular
                               params_.GetBasePolyDurations(),
                               ee_motion,
                               ee_force,
                               contact_schedule,
                               params_.IsOptimizeTimings());
  return vars;
}

std::vector<NodesVariables::Ptr>
NlpFormulationPlanar::MakeBaseVariables () const
{
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_.GetBasePolyDurations().size() + 1;

  auto spline_lin = std::make_shared<NodesVariablesPartial>(n_nodes, k3D, DimSet({X,Y}), id::base_lin_nodes);

  // double x = final_base_.lin.p().x();
  // double y = final_base_.lin.p().y();
  // double z = initial_base_.lin.p().z();
  // Vector3d final_pos = ;

  spline_lin->SetByLinearInterpolationAllNodes(initial_base_.lin.p(), final_base_.lin.p(), params_.GetTotalTime());
  spline_lin->AddStartBound(kPos, {X,Y}, initial_base_.lin.p());
  spline_lin->AddStartBound(kVel, {X,Y}, initial_base_.lin.v());
  spline_lin->AddFinalBound(kPos, params_.bounds_final_lin_pos_,   final_base_.lin.p());
  spline_lin->AddFinalBound(kVel, params_.bounds_final_lin_vel_, final_base_.lin.v());
  vars.push_back(spline_lin);

  auto spline_ang = std::make_shared<NodesVariablesPartial>(n_nodes, k3D, DimSet({Z}), id::base_ang_nodes);
  spline_ang->SetByLinearInterpolationAllNodes(initial_base_.ang.p(), final_base_.ang.p(), params_.GetTotalTime());
  spline_ang->AddStartBound(kPos, {Z}, initial_base_.ang.p());
  spline_ang->AddStartBound(kVel, {Z}, initial_base_.ang.v());
  spline_ang->AddFinalBound(kPos, params_.bounds_final_ang_pos_, final_base_.ang.p());
  spline_ang->AddFinalBound(kVel, params_.bounds_final_ang_vel_, final_base_.ang.v());
  vars.push_back(spline_ang);

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulationPlanar::MakeEndeffectorVariables () const
{
  std::vector<NodesVariablesPhaseBased::Ptr> vars;
   
  // coordinate transform to calculate nominal EE positions 
  double yaw = final_base_.ang.p().z();
  Eigen::Vector3d euler(0.0, 0.0, yaw);
  Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);

  // Endeffector Motions
  double T = params_.GetTotalTime();
  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesPlanarEEMotion>(
                                              params_.GetPhaseCount(ee),
                                              params_.ee_in_contact_at_start_.at(ee),
											  params_.min_swing_height_,
                                              id::EEMotionNodes(ee),
                                              params_.ee_polynomials_per_swing_phase_);

    // set final position
    Vector3d final_ee_W;
    if (final_ee_W_.size() == 0) {
      // use nominal final position
      final_ee_W = final_base_.lin.p() + w_R_b*model_.kinematic_model_->GetNominalStanceInBase(ee);
    }
    else {
      // use user-provided final position
      final_ee_W = final_ee_W_.at(ee);
	}
    // check if leg is in contact at finish
    bool ee_in_contact_at_finish = params_.ee_in_contact_at_start_.at(ee);
    if (params_.ee_phase_durations_.at(ee).size() % 2 == 0) ee_in_contact_at_finish = ! ee_in_contact_at_finish;
    if (ee_in_contact_at_finish) {
      // assign z to terrain level
      final_ee_W.z() = 0.0;
    }
    // initialize towards final footholds
    nodes->SetByLinearInterpolation(initial_ee_W_.at(ee), final_ee_W, T);
	// add start bounds
    nodes->AddStartBound(kPos, {X,Y}, initial_ee_W_.at(ee));
    nodes->AddStartBound(kVel, {X,Y}, Vector3d::Zero());
	// add final bounds only if they are provided in Parameters structure
    if (params_.ee_bounds_final_lin_pos_.size() != 0) {
		nodes->AddFinalBound(kPos, params_.ee_bounds_final_lin_pos_.at(ee), final_ee_W);
	}
    if (params_.ee_bounds_final_lin_vel_.size() != 0) {
		nodes->AddFinalBound(kVel, params_.ee_bounds_final_lin_vel_.at(ee), Vector3d::Zero());
	}

    vars.push_back(nodes);
  }

  return vars;
}

std::vector<NodesVariablesTimeBased::Ptr>
NlpFormulationPlanar::MakeForceVariables () const
{
  std::vector<NodesVariablesTimeBased::Ptr> vars;

  double T = params_.GetTotalTime();
  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEForceTimeBased>(
                                              params_.ee_phase_durations_.at(ee),
                                              params_.ee_in_contact_at_start_.at(ee),
                                              id::EEForceNodes(ee),
                                              params_.dt_constraint_dynamic_,
											  k1D);

    // initialize with mass of robot distributed equally on all legs
    double m = model_.dynamic_model_->m();
    double g = model_.dynamic_model_->g();

    Vector3d f_stance; //TODO optimization?
	f_stance.x() = m*g/params_.GetEECount();
    nodes->SetByLinearInterpolation(f_stance, f_stance, T); // stay constant
    vars.push_back(nodes);
  }

  return vars;
}

std::vector<PhaseDurations::Ptr>
NlpFormulationPlanar::MakeContactScheduleVariables () const
{
  std::vector<PhaseDurations::Ptr> vars;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto var = std::make_shared<PhaseDurations>(ee,
                                                params_.ee_phase_durations_.at(ee),
                                                params_.ee_in_contact_at_start_.at(ee),
                                                params_.bound_phase_duration_.first,
                                                params_.bound_phase_duration_.second);
    vars.push_back(var);
  }

  return vars;
}


NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::GetConstraints(const SplineHolder& spline_holder) const
{
  ContraintPtrVec constraints;
  for (auto name : params_.constraints_)
    for (auto c : GetConstraint(name, spline_holder))
      constraints.push_back(c);

  return constraints;
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::GetConstraint (Parameters::ConstraintName name,
                           const SplineHolder& s) const
{
  switch (name) {
    case Parameters::Dynamic:                 return MakeDynamicConstraint(s);
    case Parameters::EndeffectorRom:          return MakeRangeOfMotionConstraint(s);
    case Parameters::EndeffectorSphericalRom: return MakeRangeOfMotionSphericalConstraint(s);
    case Parameters::Force:                   return MakeForceConstraint();
    case Parameters::Swing:                   return MakeSwingConstraint();
    default: throw std::runtime_error("constraint not defined!");
  }
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::MakeDynamicConstraint(const SplineHolder& s) const
{
  auto constraint = std::make_shared<DynamicZMPPlanarConstraint>(model_.dynamic_model_,
                                                        params_.GetTotalTime(),
                                                        params_.dt_constraint_dynamic_,
                                                        s);
  return {constraint};
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::MakeRangeOfMotionConstraint (const SplineHolder& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    ContraintPtr rom;

	KinematicModel::AlignedBox3d box = model_.kinematic_model_->GetBoundingBox(ee);
	Sphere3d sphere = model_.kinematic_model_->GetBoundingSphere(ee);

    // check if bounding box is inside bounding sphere
    if ( sphere.contains(box.min()) && sphere.contains(box.max()) ) {
	  // use simplified constraint without sphere
      rom = std::make_shared<RangeOfMotionPlanarConstraint>(box, params_.GetTotalTime(), params_.dt_constraint_range_of_motion_, ee, s);
    }
    else {
	  // use constraint without sphere
      rom = std::make_shared<RangeOfMotionSphericalPlanarConstraint>(box, sphere, params_.GetTotalTime(), params_.dt_constraint_range_of_motion_, ee, s);
    }
    c.push_back(rom);
  }

  return c;
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::MakeRangeOfMotionSphericalConstraint (const SplineHolder& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto rom = std::make_shared<RangeOfMotionSphericalPlanarConstraint>(model_.kinematic_model_->GetBoundingBox(ee),
                                                         model_.kinematic_model_->GetBoundingSphere(ee),
                                                         params_.GetTotalTime(),
                                                         params_.dt_constraint_range_of_motion_,
                                                         ee,
                                                         s);
    c.push_back(rom);
  }

  return c;
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::MakeForceConstraint () const
{
  ContraintPtrVec constraints;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto c = std::make_shared<ForcePlanarConstraint>(params_.force_limit_in_normal_direction_,
                                               ee);
    constraints.push_back(c);
  }

  return constraints;
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::MakeSwingConstraint () const
{
  ContraintPtrVec constraints;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto swing = std::make_shared<SwingConstraint>(id::EEMotionNodes(ee));
    constraints.push_back(swing);
  }

  return constraints;
}

NlpFormulationPlanar::ContraintPtrVec
NlpFormulationPlanar::GetCosts(const SplineHolder& s) const
{
  ContraintPtrVec costs;
  for (const auto& pair : params_.costs_)
    for (auto c : GetCost(pair.first, pair.second, s))
      costs.push_back(c);

  return costs;
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::GetCost(const Parameters::CostName& name, double weight, const SplineHolder& s) const
{
  switch (name) {
    case Parameters::ForcesCostID:   return MakeForcesCost(weight);
    case Parameters::EEMotionCostID: return MakeEEMotionCost(weight);
    case Parameters::BaseLinMotionCostID: return MakeBaseAngMotionCost(weight);
    case Parameters::BaseAngMotionCostID: return MakeBaseLinMotionCost(weight);
    case Parameters::EEAccCostID: return MakeEEAccCost(weight, s);
    case Parameters::BaseLinAccCostID: return CostPtrVec(1, std::make_shared<AccelSquareNormCost>(id::base_lin_nodes, weight, s.base_linear_) );
    case Parameters::BaseAngAccCostID: return CostPtrVec(1, std::make_shared<AccelSquareNormCost>(id::base_ang_nodes, weight, s.base_angular_) );
	default: throw std::runtime_error("cost not defined!");
  }
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::MakeForcesCost(double weight) const
{
  CostPtrVec cost;

  for (int ee=0; ee<params_.GetEECount(); ee++)
    cost.push_back(std::make_shared<NodeCost>(id::EEForceNodes(ee), kPos, X, weight));

  return cost;
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::MakeEEMotionCost(double weight) const
{
  CostPtrVec cost;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, X, weight));
    cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, Y, weight));
  }

  return cost;
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::MakeBaseLinMotionCost(double weight) const
{
  CostPtrVec cost;

  for(int dim = X; dim <= Z; dim++) {
    cost.push_back(std::make_shared<NodeCost>(id::base_lin_nodes, kVel, dim, weight));
  }

  return cost;
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::MakeBaseAngMotionCost(double weight) const
{
  CostPtrVec cost;

  for(int dim = X; dim <= Z; dim++) {
    cost.push_back(std::make_shared<NodeCost>(id::base_ang_nodes, kVel, dim, weight));
  }

  return cost;
}

NlpFormulationPlanar::CostPtrVec
NlpFormulationPlanar::MakeEEAccCost(double weight, const SplineHolder& s) const
{
  CostPtrVec cost;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    cost.push_back(std::make_shared<AccelSquareNormCost>(id::EEMotionNodes(ee), weight, s.ee_motion_[ee]));
  }

  return cost;
}


} /* namespace towr */
