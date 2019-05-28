#include <towr/models/general_kinematic_model.h>

namespace towr {

GeneralKinematicModel::GeneralKinematicModel(int n_ee)
{
	if (n_ee > max_number_of_legs_) throw std::invalid_argument("Too many legs");
	nominal_stance_.resize(n_ee);
	bounding_box_.resize(n_ee);
	bounding_sphere_.resize(n_ee);
}

void GeneralKinematicModel::configureEndEffector(EE ee, const Vector3d& nominal_stance, const AlignedBox3d& bounding_box, const Sphere3d& bounding_sphere)
{
	int n_ee = nominal_stance_.size();
	if (ee >= n_ee) throw std::out_of_range("invalid end effector number");
	// check parameters sanity
	if (!bounding_box.contains(nominal_stance) || !bounding_sphere.contains(nominal_stance)) {
		throw std::invalid_argument("nominal pose outside range restrictions");
	}
	// assign end effector parameters
	nominal_stance_[ee] = nominal_stance;
	bounding_box_[ee] = bounding_box;
	bounding_sphere_[ee] = bounding_sphere;
	configured_legs_[ee] = true;
}

} /* namespace towr */


