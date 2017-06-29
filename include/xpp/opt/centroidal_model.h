/**
 @file    centroidal_model.h
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    May 19, 2017
 @brief   Brief description
 */

#ifndef XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_
#define XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_

#include "dynamic_model.h"

namespace xpp {
namespace opt {

/**
 * @brief Centroidal Dynamics for the 6-DoF Base to model the system.
 */
class CentroidalModel : public DynamicModel {
public:
  CentroidalModel ();
  virtual ~CentroidalModel ();

  virtual BaseAcc GetBaseAcceleration() const override;

  virtual Jacobian GetJacobianOfAccWrtBaseLin(const BaseLin&,
                                              double t_global) const override;
  virtual Jacobian GetJacobianOfAccWrtBaseAng(const BaseAng&,
                                            double t_global) const override;
  virtual Jacobian GetJacobianofAccWrtForce(const EndeffectorsForce&,
                                            double t_global,
                                            EndeffectorID) const override;
  virtual Jacobian GetJacobianofAccWrtEEPos(const EndeffectorsMotion&,
                                             double t_global,
                                             EndeffectorID) const override;
private:
  double m_;          /// mass of robot
  Eigen::Matrix3d I_; /// inertia tensor of robot
  Eigen::SparseMatrix<double, Eigen::RowMajor> I_inv_;


  static Jacobian
  BuildCrossProductMatrix(const Vector3d& in)
  {
    Jacobian out(3,3);

    out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
    out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
    out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

    return out;
  }

  static Eigen::Matrix3d buildInertiaTensor(
          double Ixx, double Iyy, double Izz,
          double Ixy, double Ixz, double Iyz)
  {
    Eigen::Matrix3d I;
    I <<  Ixx, -Ixy, -Ixz,
         -Ixy,  Iyy, -Iyz,
         -Ixz, -Iyz,  Izz;
    return I;
  }
};

} /* namespace opt */
} /* namespace xpp */

#endif /* XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_ */