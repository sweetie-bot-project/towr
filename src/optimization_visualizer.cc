/**
 @file    optimization_visualizer.cc
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    May 31, 2016
 @brief   Brief description
 */

#include <xpp/ros/optimization_visualizer.h>
#include <xpp/ros/ros_helpers.h>

namespace xpp {
namespace ros {

OptimizationVisualizer::OptimizationVisualizer ()
{
  InitRos();
}

OptimizationVisualizer::OptimizationVisualizer (OptimizationVariables& subject)
    : subject_(&subject)
{
  RegisterWithSubject(subject);
  InitRos();
}

void
OptimizationVisualizer::RegisterWithSubject (OptimizationVariables& subject)
{
  subject_ = &subject;
  subject_->RegisterObserver(this);
}

void
OptimizationVisualizer::SetInterpreter (const Interpreter& interpreter)
{
  interpreter_ = interpreter;
}

void
OptimizationVisualizer::Update ()
{
  splines_   = interpreter_.GetSplines(subject_->GetSplineCoefficients());
  footholds_ = interpreter_.GetFootholds(subject_->GetFootholdsStd());
}

void
OptimizationVisualizer::PublishMsg ()
{
  double walking_height = RosHelpers::GetDoubleFromServer("/xpp/robot_height");

  visualization_msgs::MarkerArray msg = msg_builder_.BuildMsg(splines_, footholds_, walking_height);
  ros_publisher_.publish(msg);
}

void
OptimizationVisualizer::InitRos ()
{
  ::ros::NodeHandle n;
  ros_publisher_ = n.advertise<visualization_msgs::MarkerArray>("optimization_variables", 1);
}

} /* namespace ros */
} /* namespace xpp */
