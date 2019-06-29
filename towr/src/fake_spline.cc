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

#include <towr/variables/fake_spline.h>

namespace towr {

FakeSpline::FakeSpline(NodesVariablesTimeBased * node_variables)
    : NodesObserver(node_variables)
{
	node_values_time_based_ = node_variables;
} 

void
FakeSpline::UpdateNodes ()
{
  // do nothing: all information is stored in nodes.
}

int
FakeSpline::GetNodeVariablesCount() const
{
  return node_values_->GetRows();
}

const State FakeSpline::GetPoint(double t) const
{
  int node_id = node_values_time_based_->GetNodeIdExact(t);
  const std::vector<Node>& nodes = node_values_->GetNodes();
  return nodes.at(node_id);
}

FakeSpline::Jacobian
FakeSpline::GetJacobianWrtNodes (double t_global, Dx dxdt) const
{
  assert(dxdt == kPos);
  int node_id = node_values_time_based_->GetNodeIdExact(t_global);
  Jacobian jac = Jacobian(node_values_->GetDim(), node_values_->GetRows());

  for(int dim = 0; dim < k3D; dim++) {
    int idx = node_values_->GetOptIndex(NodesVariables::NodeValueInfo(node_id, kPos, dim));
	if (idx != NodesVariables::NodeValueNotOptimized) {
	  jac.coeffRef(dim, idx) = 1.0;
	}
  }

  return jac;
}

FakeSpline::Jacobian
FakeSpline::GetJacobianWrtNodes (int id, double t_local, Dx dxdt) const
{
  assert(dxdt == kPos);

  // calculate and check node_id
  const double eps = 1e-10; // double precision
  int node_id;
  if (std::abs(t_local) < eps) {
    node_id = id;
  }
  else {
	node_id = id+1;
	double dT = node_values_time_based_->GetDiscretizationTime();
	assert(std::abs(t_local - dT) < eps);
  }

  Jacobian jac = Jacobian(node_values_->GetDim(), node_values_->GetRows());

  for(int dim = 0; dim < k3D; dim++) {
    int idx = node_values_->GetOptIndex(NodesVariables::NodeValueInfo(node_id, kPos, dim));
	if (idx != NodesVariables::NodeValueNotOptimized) {
	  jac.coeffRef(dim, idx) = 1.0;
	}
  }

  return jac;
}

} /* namespace towr */
