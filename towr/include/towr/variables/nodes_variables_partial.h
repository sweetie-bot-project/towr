/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. Partial rights reserved.

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

#ifndef TOWR_TOWR_INCLUDE_TOWR_VARIABLES_PARTIAL_BASE_NODES_H_
#define TOWR_TOWR_INCLUDE_TOWR_VARIABLES_PARTIAL_BASE_NODES_H_

#include "nodes_variables.h"

namespace towr {

/**
 * @brief Node variables used to construct the base motion spline. For some dimentions node values are fixed.
 *
 * Every node is optimized overlike NodeVariablesAll, but some dimentions are fixed and not being optimized.
 * This is useful if planar or another some kind of restricted motion is modelled.
 *
 * @ingroup Variables
 */
class NodesVariablesPartial : public NodesVariables {
public:
  /**
   * @param n_nodes  Number of nodes to construct the spline.
   * @param n_dim    Number of dimensions of each node.
   * @param optim_dims    Set of dimentions being optimized.
   * @param variable_id  Name of this variables set in the optimization.
   */
  NodesVariablesPartial (int n_nodes, int n_dim, DimSet optim_dims, std::string variable_id);
  virtual ~NodesVariablesPartial () = default;

  std::vector<NodeValueInfo> GetNodeValuesInfo(int idx) const override;
  int GetOptIndex(const NodeValueInfo& nvi) const override;

protected:
  DimSet optim_dims_;
};

} /* namespace towr */

#endif /* TOWR_TOWR_INCLUDE_TOWR_VARIABLES_PARTIAL_BASE_NODES_H_ */