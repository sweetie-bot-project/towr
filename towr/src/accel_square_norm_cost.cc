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

#include <towr/costs/accel_square_norm_cost.h>
#include <iostream>

namespace towr {

AccelSquareNormCost::AccelSquareNormCost (const std::string& nodes_id, double weight, const NodeSpline::Ptr spline)
    : CostTerm(nodes_id+"-accel_norm")
{
  node_id_ = nodes_id;
  weight_ = weight;
  spline_ = spline;
}

double
AccelSquareNormCost::GetCost () const
{
  return weight_ * spline_->GetAccSquareNormValue().sum();
}

void
AccelSquareNormCost::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == node_id_) {
	Jacobian jac_full = weight_ * spline_->GetAccSquareNormJacobianWrtNodes();
	for(int k = 0; k < jac_full.rows(); k++) jac += jac_full.row(k);

	//std::cout << "\nVARSET: " << var_set << "\ngrad3:\n" << jac_full << "\ngrad3:\n" << jac << "\n";
  }
}

} /* namespace towr */

