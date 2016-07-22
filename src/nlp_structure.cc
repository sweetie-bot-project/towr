/**
 @file    nlp_structure.cc
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    Jun 8, 2016
 @brief   Defines the class functions in NlpStructure
 */

#include <xpp/zmp/nlp_structure.h>

namespace xpp {
namespace zmp {

class VariableSet {
public:
  typedef Eigen::VectorXd VectorXd;
  typedef AConstraint::VecBound VecBound;

  VariableSet(int n_variables, const std::string& name);
  virtual ~VariableSet();

  VectorXd GetVariables() const;
  VecBound GetBounds() const;
  std::string GetName() const;

  void SetVariables(const VectorXd& x);
private:
  VectorXd x_;
  VecBound bounds_;
  std::string name_;
};

VariableSet::VariableSet (int n_variables, const std::string& name)
{
  x_ = Eigen::VectorXd::Zero(n_variables);
  bounds_.assign(n_variables, AConstraint::kNoBound_);
  name_ = name;
}

VariableSet::~VariableSet ()
{
}

std::string
VariableSet::GetName () const
{
  return name_;
}

VariableSet::VectorXd
VariableSet::GetVariables () const
{
  return x_;
}

VariableSet::VecBound
VariableSet::GetBounds () const
{
  return bounds_;
}

void
VariableSet::SetVariables (const VectorXd& x)
{
  x_ = x;
}
///////////////////////////////////////////////////////////////////////////////

NlpStructure::NlpStructure()
{
  n_variables_ = 0;
}

NlpStructure::~NlpStructure ()
{
}

void
NlpStructure::AddVariableSet (std::string name, int n_variables)
{
  variable_sets_.push_back(VariableSetPtr(new VariableSet(n_variables, name)));
  n_variables_ += n_variables;
}

void
NlpStructure::Reset ()
{
  variable_sets_.clear();
  n_variables_ = 0;
}

int
NlpStructure::GetOptimizationVariableCount() const
{
  return n_variables_;
}

NlpStructure::VectorXd
NlpStructure::GetAllOptimizationVariables () const
{
  Eigen::VectorXd x(GetOptimizationVariableCount());
  int c = 0;
  for (const auto& set : variable_sets_) {
    const VectorXd& var = set->GetVariables();
    x.middleRows(c, var.rows()) = var;
    c += var.rows();
  }

  return x;
}

NlpStructure::VecBound
NlpStructure::GetAllBounds () const
{
  VecBound bounds_;
  for (const auto& set : variable_sets_) {
    const VecBound& b = set->GetBounds();
    bounds_.insert(std::end(bounds_), std::begin(b), std::end(b));
  }

  return bounds_;
}

void
NlpStructure::SetAllVariables(const VectorXd& x_all)
{
  int c = 0;
  for (const auto& set : variable_sets_) {
    int n_var = set->GetVariables().rows();
    set->SetVariables(x_all.middleRows(c,n_var));
    c += n_var;
  }
}

void
NlpStructure::SetVariables (std::string set_name, const VectorXd& values)
{
  for (const auto& s : variable_sets_)
    if (s->GetName() == set_name) {
     s->SetVariables(values);
     return;
    }

  assert(false); // name not present in set
}

NlpStructure::VectorXd
NlpStructure::GetVariables (std::string set_name) const
{
  for (const auto& s : variable_sets_)
    if (s->GetName() == set_name)
     return s->GetVariables();

  assert(false); // name not present in set
}

} // namespace zmp
} // namespace xpp