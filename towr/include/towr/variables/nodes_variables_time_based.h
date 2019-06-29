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

#ifndef TOWR_TOWR_INCLUDE_TOWR_VARIABLES_TIME_BASED_H_
#define TOWR_TOWR_INCLUDE_TOWR_VARIABLES_TIME_BASED_H_

#include "nodes_variables.h"

namespace towr {

class NodesVariablesTimeBased : public NodesVariables {
public:
  using Ptr         = std::shared_ptr<NodesVariablesTimeBased>;
  using NodeIds     = std::vector<int>;
  using OptIndex    = std::vector<NodeValueInfo>; 
  using VecTimes    = std::vector<double>;

protected:
  /**
   * @param variable_name  Name of this variables set in the optimization.
   * @param dT discretization time. 
   * @param T Total duration.
   */
  NodesVariablesTimeBased (const std::string& variable_name, double dT, double T) 
	  : NodesVariables(variable_name), dT_(dT), TotalTime_(T) {}

public:
  virtual ~NodesVariablesTimeBased () = default;

  VecTimes GetTimeIntervalDurations() const;
  int GetNodeIdExact(double t_global) const;
  /**
   * @brief Return time instance, associated with node.
   * @param node_id Node index.
   * @return Time instance assocated with node.
   **/
  double GetNodeTime(int node_id) const;

  double GetTotalDuraion() const { return TotalTime_; };
  double GetDiscretizationTime() const { return dT_; };

protected:
  double TotalTime_ /**< Total time duration. */;
  double dT_; /**< Time discretization interval duration */
};

/**
 * @brief Node variables for representing reaction forses. It designed to match DynamicConstraint with less overheat.
 *
 * One node match to each constraint in TimeBasedConstraint. If constraints is equalities this significantly
 * reduces computation cost of matching to splines via constrints.  Only nodes' positions are optimized. 
 * During swing phases forces is asumed equal to zero. 
 *
 * @ingroup Variables
 */
class NodesVariablesEEForceTimeBased : public NodesVariablesTimeBased {
public:
  using Ptr         = std::shared_ptr<NodesVariablesEEForceTimeBased>;

  struct NodeInfo {
    bool is_constant_; /**< if node is not optimized */
	int phase_; /**< node's phase */
	int opt_index_base_; /** Start of block of opt iducies assotiated with node. */

	NodeInfo() = default;
	NodeInfo(bool is_constant, int phase, int opt_index_base) 
		: is_constant_(is_constant), phase_(phase), opt_index_base_(opt_index_base) {} 
  };

public:
  /**
   * @param phase_durations Durations of motion phases for end effector.
   * @param first_phase_contact True if end effector is in contact at start.
   * @param variable_id  Name of this variables set in the optimization.
   * @param dT discretization time. Should exactly match discretization time of coresponding TimeBasedConstraint.
   */
  NodesVariablesEEForceTimeBased(const VecTimes& phase_durations, bool first_phase_contact, std::string variable_id, double dT);
  virtual ~NodesVariablesEEForceTimeBased () = default;

  /**
   * @brief The indices of those nodes that don't belong to a constant phase.
   *
   * For forces nodes these are the stance phases (can produce force), and for
   * feet these are the swing phases (can move end-effector).
   */
  NodeIds GetIndicesOfNonConstantNodes() const;
  int GetPhase(int node_id) const;

  int GetOptIndex(const NodeValueInfo& nvi) const override;
  std::vector<NodeValueInfo> GetNodeValuesInfo(int idx) const override;

private:
  std::vector<NodeInfo> nodes_info_;  // node_id -> information about node
  OptIndex index_to_node_value_info_; // opt index -> information about node
};

} /* namespace towr */

#endif /* TOWR_TOWR_INCLUDE_TOWR_VARIABLES_TIME_BASED_H_ */
