/*!
 * \file   zmp_runner.cc
 * \author Alexander Winkler
 * \date   Oct 4, 2014
 * \brief  SL Task that executes a walking gait given a arbitrary sequence and
 *         position of footholds.
 *
 *         It uses the dynamic locomotion library \c xpp to find an optimal and
 *         dynamicaly stable body trajectory and executes this trjaectory using
 *         inverse dynamics coupled with a virtual model controller.
 */


#include <xpp/exe/walking_controller.h>

#include <InverseKinematics.h>
#include <xpp/ros/ros_helpers.h>
#include <xpp/exe/robot_interface.h>

#include <kindr/rotations/RotationEigen.hpp>


namespace xpp {
namespace exe {

// only for logging
static Eigen::Vector3d log_rpy_curr;
static Eigen::Vector3d log_rpy_des;
static int log_swingleg_id;
static xpp::hyq::VirtualModel::Accelerations log_base_acc_des_fb;
static xpp::hyq::VirtualModel::Accelerations log_base_acc_des_ff;


WalkingController::WalkingController()
     :jsim_(inertia_properties_, force_transforms_)
{
}

WalkingController::~WalkingController()
{
}

void
WalkingController::GetReadyHook() {
  using namespace xpp::ros;
  // setup ROS
  int argc = 0;
  char** argv;
  ::ros::init(argc, argv, "sl_ros_subscriber_node");
  std::shared_ptr<::ros::NodeHandle> nh(new ::ros::NodeHandle);
  current_info_pub_ = nh->advertise<ReqInfoMsg>("required_info_nlp", 1);
  opt_params_sub_ = nh->subscribe("optimized_parameters_nlp", 1, &WalkingController::OptParamsCallback, this);


  b_r_geomtocog = inertia_properties_.getCOM_trunk();
  //  b_r_geomtocog(Y) += 0.020; // because of hoses

  use_virtual_model_ = RosHelpers::GetBoolFromServer("/exec/use_virtual_model");
  vm_.SetGains(500,   1200,  800, // position gains x,y,z -> drift    | 200, 600, 1000,
               200,   200,  400, // velocity gains x,y,z             | 200, 200, 400,
               700,   700, 1000, // position gains roll, pitch, yaw  | 1000, 1000, 1000
               100,   100,  200);// velocity gains roll, pitch, yaw  | 100, 100, 100);

  ffspline_duration_ = RosHelpers::GetDoubleFromServer("/exec/ff_spline_duration");

  double lift_height = RosHelpers::GetDoubleFromServer("/exec/lift_height");
  double outward_swing = RosHelpers::GetDoubleFromServer("/exec/outward_swing_distance");
  spliner_.SetParams(0.5, lift_height, outward_swing);

  robot_height_     = RosHelpers::GetDoubleFromServer("/xpp/robot_height");

  // to determine when to start reoptimization
  t_stance_initial_ = RosHelpers::GetDoubleFromServer("/xpp/stance_time_initial");
  t_swing_          = RosHelpers::GetDoubleFromServer("/xpp/swing_time");

  states_map_ = WalkingControllerState::BuildStates();
  current_state_ = WalkingControllerState::kFirstPlanning;

 first_time_sending_commands_ = true;
}

bool
WalkingController::RunHook()
{
  // delegates the actual execution to the current state
  states_map_.at(current_state_)->Run(this);
  return true;
}

void
WalkingController::SetState(WalkingControllerState::State state) {
  current_state_ = state;
}

void
WalkingController::OptParamsCallback(const OptimizedParametersMsg& msg)
{
  opt_splines_   = xpp::ros::RosHelpers::RosToXpp(msg.splines);
  opt_footholds_ = xpp::ros::RosHelpers::RosToXpp(msg.footholds);

  ROS_INFO_STREAM("updated splines [size=" << opt_splines_.size() << "] and footholds [size=" << opt_footholds_.size() << "]");
}

void WalkingController::PublishCurrentState()
{
  //    AddVarForLogging();

  ReqInfoMsg msg;
  msg.curr_state    = xpp::ros::RosHelpers::XppToRos(P_curr_.base_.pos);
  msg.curr_stance   = xpp::ros::RosHelpers::XppToRos(P_curr_.GetStanceLegs());
  msg.curr_swingleg = xpp::hyq::NO_SWING_LEG;

  // send out the message
  current_info_pub_.publish(msg);

  P_des_ = P_curr_;
}

void WalkingController::BuildPlan()
{
  ffspliner_timer_ = ffspline_duration_;
  reoptimize_before_finish_ = true;

  ::ros::spinOnce(); // process callbacks (get the optimized values).

  // start from desired state so there is no jump in desired
  spliner_.Init(P_des_, opt_splines_, opt_footholds_, robot_height_);


  kOptTimeReq_ = 0.45; // allow >150ms for ROS communication. how long before end of spline to start optimizing
  switch_node_ = spliner_.GetNode(1);

  t_switch_ = switch_node_.T;
  Controller::ResetTime();
}


void WalkingController::ExecuteLoop()
{
  using namespace xpp::utils;
  using namespace xpp::zmp;
  /** 1. motion plan generation */
  /** CURRENT state of robot through joint encoder readings and state estimation **/
  JointState q = robot_->GetJointPosition();
  JointState qd = robot_->GetJointVelocity();

  P_curr_.swingleg_ = P_des_.swingleg_;
  EstimateCurrPose(); // through sensors and state estimation
  jsim_.update(q);

  double time_left = t_switch_ - Time();
  if (time_left <= kOptTimeReq_ && reoptimize_before_finish_)
  {
    Point3d predicted_state = GetStartStateForOptimization();

    // todo this is the one that doesn't let me change nodes at arbitrary times (stance different from 4 legs)
    // use this somehow:
    //    LegDataMap<Point3d> predicted_feet;
    //    LegDataMap<bool> predicted_swing_leg;
    //    spliner_.FillCurrFeet(t_switch_, predicted_feet, predicted_swing_leg);
    // or VecFoothold predicted_stance = spliner_.GetGoalNode(Time()).state_.GetStanceLegs();
    VecFoothold predicted_stance = switch_node_.state_.FeetToFootholds().ToVector();


    ROS_INFO_STREAM("time: " << Time() << "\npredicted_start_state_:\n" << predicted_state);

    std::cout << "predicted stance: \n";
    for (const hyq::Foothold& f : predicted_stance)
      std::cout << f << std::endl;

    // always just send current information, no logic/commands/etc
    ReqInfoMsg msg;
    msg.curr_state = xpp::ros::RosHelpers::XppToRos(predicted_state);
    msg.curr_stance = xpp::ros::RosHelpers::XppToRos(predicted_stance);
    msg.curr_swingleg = P_curr_.SwinglegID();

    current_info_pub_.publish(msg);

    reoptimize_before_finish_ = false;
  }




  /** @brief DESIRED state given by splined plan and zmp optimizer **/
  P_des_.base_.pos = spliner_.GetCurrPosition(Time());
  P_des_.base_.ori = spliner_.GetCurrOrientation(Time());
  spliner_.FillCurrFeet(Time(), P_des_.feet_, P_des_.swingleg_);
  // logging
  log_base_acc_des_ff.segment<3>(LX) = P_des_.base_.pos.a; // logging only
  ROS_INFO_STREAM_THROTTLE(robot_->GetControlLoopInterval(), "time: " << Time());
//  ROS_INFO_STREAM_THROTTLE(dt_, "\nP_des_:\n" << P_des_);
  Orientation::QuaternionToRPY(P_des_.base_.ori.q, log_rpy_des);
  log_swingleg_id = P_des_.SwinglegID();




  /** TRACKING of desired body motion */
  VirtualModel::Accelerations P_base_acc_des;
  P_base_acc_des.segment(AX, 3) =  P_des_.base_.ori.a;
  P_base_acc_des.segment(LX, 3) =  P_des_.base_.pos.a;
  log_base_acc_des_ff = P_base_acc_des; // logging only

  if (use_virtual_model_) {
    // add feed-back accelerations if deviation from desired pose
    VirtualModel::Accelerations i_base_acc_des_fb;

    vm_.CalcFeedbackAcc(P_curr_.base_, P_des_.base_, jsim_.getWholeBodyInertia(),
                        i_base_acc_des_fb);
    log_base_acc_des_fb = i_base_acc_des_fb; // logging only

    P_base_acc_des += i_base_acc_des_fb;
  }



  // transform to desired joint state
  JointState q_des;
  HyQInverseKinematics inverseKinematics;
  Eigen::Vector3d q_des_leg;
  Eigen::Matrix3d P_R_Bdes = P_des_.base_.ori.q.normalized().toRotationMatrix();
  for (size_t ee=0; ee<xpp::hyq::_LEGS_COUNT; ee++) {

    // Transform global into local feet position
    Eigen::Vector3d P_base = P_des_.base_.pos.p;
    Eigen::Vector3d P_foot = P_des_.GetFeetPosOnly()[ee];
    Eigen::Vector3d B_foot = P_R_Bdes.transpose() * (P_foot - P_base);

    inverseKinematics.update(ee, B_foot);
    inverseKinematics.getJointPosition(q_des_leg);
    q_des.segment<3>(3*ee) = q_des_leg;
  }


  // fixme: maybe move this to robot scope
  JointState qd_des = robot_->EstimateDesiredJointVelocity(q_des, first_time_sending_commands_);
  JointState qdd_des = robot_->EstimateDesiredJointAcceleration(qd_des, first_time_sending_commands_);


  JointState uff = robot_->CalcRequiredTorques(q_des, qd_des, qdd_des,P_base_acc_des,P_curr_.swingleg_);


  SmoothTorquesAtContactChange(uff);

  if (std::fabs(uff.maxCoeff()) > 180)
    throw std::runtime_error("uff torques greater than 150");

  // write PD torque references and feed-forward torques in appropriate SL variable
  robot_->SetDesiredJointPosition(q_des);
  robot_->SetDesiredJointVelocity(qd_des);
  robot_->SetDesiredTorque(uff);
  first_time_sending_commands_ = false;
}

WalkingController::State
WalkingController::GetStartStateForOptimization(/*const double required_time*/) const
{
  // fixme initialize not with planned, but with current predicted state
  // this allows to react to pushes/disturbances
  using namespace xpp::utils;
  using namespace xpp::zmp;

////  // predict where current state is going to be after optimization is done (e.g. after kOptReq)
//  double opt_time = 0.3;
//  State curr_state = P_curr_.base_.pos;
//  ROS_INFO_STREAM("curr_state: \t" << curr_state);
//
//  Eigen::Vector3d jerk = (curr_state.a - prev_state_.a)/robot_->GetControlLoopInterval();
//
//  CoeffValues coeff;
//  coeff.x[F] = curr_state.p.x();
//  coeff.x[E] = curr_state.v.x();
//  coeff.x[D] = curr_state.a.x()/2.;
//  coeff.x[C] = jerk.x()/6.;
//
//  coeff.y[F] = curr_state.p.y();
//  coeff.y[E] = curr_state.v.y();
//  coeff.y[D] = curr_state.a.y()/2.;
//  coeff.y[C] = jerk.y()/6.;
//
//  Spline spline(coeff);
//
//  State predicted_state = curr_state; // z position, vel and acceleration
//  predicted_state.p.segment<kDim2d>(X) = spline.GetState(kPos, opt_time);
//  predicted_state.v.segment<kDim2d>(X) = spline.GetState(kVel, opt_time);
//  predicted_state.a.segment<kDim2d>(X) = spline.GetState(kAcc, opt_time);
//  ROS_INFO_STREAM("predicted: \t" << predicted_state);
//  return predicted_state;


  // this is the state where the robot is planned to be, no current feedback
  State end_des;
  Point2d end_des_xy = spliner_.GetCurrPosition(t_switch_).Get2D();
  end_des.p.segment(0,2) = end_des_xy.p; // - b_r_geomtocog.segment<2>(X)*/;
  end_des.v.segment(0,2) = end_des_xy.v;
  end_des.a.segment(0,2) = end_des_xy.a;
  return end_des;
}




void WalkingController::EstimateCurrPose()
{
  xpp::utils::Pose curr_base_state_est = robot_->GetBodyPose();
//  GetStateEst(::base_state, ::base_orient, curr_base_state_est);


  // position and orientation are expressed in inertial frame = projected frame
  P_curr_.base_.pos.p.x() = curr_base_state_est.pos.p.x();
  P_curr_.base_.pos.p.y() = curr_base_state_est.pos.p.y();
  P_curr_.base_.ori.q = curr_base_state_est.ori.q;

  // vel/acc are expressed in body frame, transfer to projected frame
  Eigen::Matrix3d P_R_B = curr_base_state_est.ori.q.normalized().toRotationMatrix();
  P_curr_.base_.pos.v = P_R_B*curr_base_state_est.pos.v;
  P_curr_.base_.pos.a = P_R_B*curr_base_state_est.pos.a;
  P_curr_.base_.ori.v = P_R_B*curr_base_state_est.ori.v;
  P_curr_.base_.ori.a = P_R_B*curr_base_state_est.ori.a;
  P_curr_.base_.pos.a.z() -= iit::rbd::g; // exclude gravity from base acceleration


  // state estimator drifts, so estimate body position by joint angles.
  LegDataMap<Vector3d> curr_feet_local = robot_->GetFeetPositions();
  double z_sum = 0.0;
  int stance_leg_count = 0;
  for (LegID leg : hyq::LegIDArray)
  {
    if (!P_curr_.swingleg_[leg]) // only use stance legs
    {
      z_sum += curr_feet_local[leg].z();
      ++stance_leg_count;
    }
  }
  double b_tz_basetofeet = z_sum / stance_leg_count;
  Eigen::Vector3d B_t_basetofeet = {0.0, 0.0, b_tz_basetofeet};
  Eigen::Vector3d P_t_basetofeet = P_R_B * B_t_basetofeet;
  double base_new = - P_t_basetofeet.z() + P_des_.GetZAvg();


  // todo this still estimates high current z-velocities, which aren't really there
  if (first_time_sending_commands_) {
    P_curr_.base_.pos.p.z() = base_new;
  } else {
    // inifinite impulse response filter
    P_curr_.base_.pos.p.z() = 0.97*P_curr_.base_.pos.p.z() + 0.03*base_new;
  }

  // foot position must be in the same frame as the base
  for (LegID leg : hyq::LegIDArray)
    P_curr_.feet_[leg].p = TransformBaseToProjectedFrame(curr_feet_local[leg], P_curr_.base_);


  // safety check
  Orientation::QuaternionToRPY(P_curr_.base_.ori.q, log_rpy_curr);
  if (first_time_sending_commands_) {
    if(    (std::abs(log_rpy_curr.x()) > 7./180. * M_PI)
        || (std::abs(log_rpy_curr.y()) > 7./180. * M_PI)
        || (P_curr_.base_.pos.p.z() < 0.4)
    ) {
      ROS_WARN_STREAM("i_curr_: \n" << P_curr_);
      throw std::runtime_error("initial base state is not in projected frame");
    }
  }

  // logging
  ROS_DEBUG_STREAM_THROTTLE(robot_->GetControlLoopInterval(), "time: " << Time() << "\nP_curr_:\n" << P_curr_);
}

bool
WalkingController::TimeExceeded () const
{
  double t_max = t_switch_ - robot_->GetControlLoopInterval();
  return Time() >= t_max; /*|| Time() >= spliner_.GetTotalTime()-dt_*/
}

Eigen::Vector3d
WalkingController::TransformBaseToProjectedFrame(const Eigen::Vector3d& B_r_btox,
                                         const xpp::utils::Pose& P_base_BtoP) const
{
  Eigen::Matrix3d P_R_B    = P_base_BtoP.ori.q.normalized().toRotationMatrix();
  Eigen::Vector3d P_r_ptob = P_base_BtoP.pos.p;

  Eigen::Vector3d P_r_ptox = P_r_ptob + P_R_B*B_r_btox;
  return P_r_ptox;
}


void WalkingController::SmoothTorquesAtContactChange(JointState& uff)
{

  static bool ffsplining_ = false; // the first control loop ever we want no splining, full torques!

  if (ffspline_duration_ > 0.0) {
    // check if contacts have changed during this task loop
    for (LegID leg : hyq::LegIDArray) {
      if((P_curr_.swingleg_[leg] !=  prev_swingleg_[leg])
          /* &&  prev_swingleg_[leg]*/ ) {
        ffsplining_ = true;
        ffspliner_timer_ = ffspline_duration_;
        break;
      }
    }

    if (ffsplining_) {
      if (ffspliner_timer_ >= 0.0) {
        uff = uff_prev_ + ((ffspline_duration_ - ffspliner_timer_)/ffspline_duration_) * (uff - uff_prev_);
        ffspliner_timer_ -= robot_->GetControlLoopInterval();
      } else {
        ffsplining_ = false;
      }
    }
  }

  uff_prev_ = uff;
  prev_swingleg_ = P_curr_.swingleg_;
}


} // namespace exe
} // namespace xpp
