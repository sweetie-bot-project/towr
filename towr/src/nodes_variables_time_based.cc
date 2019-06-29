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

#include <towr/variables/nodes_variables_time_based.h>

#include <numeric> //accumulate

namespace towr {

int NodesVariablesTimeBased::GetNodeIdExact(double t_global) const
{
  // TODO rewrite more sensible?
  const double eps = 1e-10; // double precision
  unsigned int n_nodes = nodes_.size();
  unsigned int node_id = std::floor( (t_global + eps) / dT_ );
  if (node_id >= 0 && node_id < n_nodes-2) {
    assert( std::abs(node_id*dT_ - t_global) < eps );
	return node_id;
  }
  else if (node_id == n_nodes-2) {
    if (std::abs(TotalTime_ - t_global) < eps) {
      return n_nodes-1;
    }
	else {
      assert( std::abs(node_id*dT_ - t_global) < eps );
	  return node_id;
    }
  } 
  else if (node_id == n_nodes-1) {
    assert(std::abs(TotalTime_ - t_global) < eps);
	return node_id;
  }
  assert(false);
  return 0;
}

double NodesVariablesTimeBased::GetNodeTime(int node_id) const
{
  if (node_id >= nodes_.size()-1) return TotalTime_;
  else return node_id * dT_;
}

NodesVariablesTimeBased::VecTimes 
NodesVariablesTimeBased::GetTimeIntervalDurations() const
{
  VecTimes intervals;
  int n_full_intervals = nodes_.size() - 2;
  intervals.reserve(n_full_intervals + 1);
  for(int k = 0; k < n_full_intervals; k++) intervals.push_back(dT_);
  intervals.push_back(TotalTime_ - dT_*n_full_intervals);
  return intervals;
}


NodesVariablesEEForceTimeBased::NodesVariablesEEForceTimeBased(const VecTimes& phase_durations, bool first_phase_contact, std::string variable_id, double dT) 
    : NodesVariablesTimeBased(variable_id, dT, std::accumulate(phase_durations.begin(), phase_durations.end(), 0.0) )
{
  // total time
  int n_nodes = std::floor(TotalTime_ / dT_)+2;

  // setup nodes list
  n_dim_ = k3D;
  nodes_ = std::vector<Node>(n_nodes, Node(n_dim_));
  nodes_info_.reserve(n_nodes);

  int n_opt_variables = 0;
  int phase_idx = 0; // index motion phase
  double phase_local_t; // time from start of phase
  bool phase_contact = first_phase_contact; // true if current phase is contact

  for (int node_id=0; node_id<nodes_.size(); ++node_id) {
    if (phase_contact) {
	  // stance node:
      // forces can be created during stance, so these nodes are optimized over.
      for (int dim=0; dim<n_dim_; ++dim) {
		// only position is optimized
        index_to_node_value_info_.emplace_back(node_id, kPos, dim); // add NodeValueInfo to index
        nodes_.at(node_id).at(kVel).setZero();
      }
	  // add metainformation about node
	  nodes_info_.emplace_back(false, phase_idx, n_opt_variables);
	  // increase number of optimized variables
	  n_opt_variables += n_dim_;
	}
    else {
	  // swing node 
      // forces can't exist during swing phase, so no need to be optimized
      // -> all node values simply set to zero.
      nodes_.at(node_id).at(kPos).setZero();
      nodes_.at(node_id).at(kVel).setZero();
	  // add metainformation about node
	  nodes_info_.emplace_back(true, phase_idx, NodeValueNotOptimized);
    }
	// shift time and calculate next phase_contact value
	phase_local_t += dT_;
	while (phase_local_t > phase_durations[phase_idx]) {
	  if (phase_idx < phase_durations.size()-1) {
		// current phase has ended: toggle phase_contact and switch to next phase
		phase_contact = !phase_contact;
		phase_local_t -= phase_durations[phase_idx];
	    phase_idx++;
	  }
	  else break;
	}
  }

  // setup otimization variables
  bounds_ = VecBound(n_opt_variables, ifopt::NoBound);
  SetRows(n_opt_variables);
}

int 
NodesVariablesEEForceTimeBased::GetPhase(int node_id) const
{
  return nodes_info_.at(node_id).phase_;
}

NodesVariablesEEForceTimeBased::NodeIds
NodesVariablesEEForceTimeBased::GetIndicesOfNonConstantNodes() const
{
  NodeIds node_ids;
  // use index to get node_id of non-constatant nodes
  for (int id = 0; id < nodes_info_.size(); id++) {
    if (!nodes_info_[id].is_constant_) node_ids.push_back(id);
  }
  return node_ids;
}

int 
NodesVariablesEEForceTimeBased::GetOptIndex(const NodeValueInfo& nvi) const 
{
  assert(nvi.deriv_ == kPos);
  const NodeInfo& ni = nodes_info_.at(nvi.id_);
  if (!ni.is_constant_) return ni.opt_index_base_ + nvi.dim_;
  else return NodeValueNotOptimized;
}

std::vector<NodesVariablesEEForceTimeBased::NodeValueInfo>
NodesVariablesEEForceTimeBased::GetNodeValuesInfo (int idx) const
{
  std::vector<NodeValueInfo> nodes_info;
  nodes_info.push_back( index_to_node_value_info_.at(idx) );
  return nodes_info;
}

} /* namespace towr */
