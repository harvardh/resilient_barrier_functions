//
//  barrier.h
//  leaderfollower
//
//  Created by Harvard Virgil Humphrey on 2019/7/15.
//  Copyright © 2019 Harvard Virgil Humphrey. All rights reserved.
//

#ifndef barrier_h
#define barrier_h
#include <thread>
#include <chrono>
#include <random>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <mav_msgs/conversions.h>
#include <mav_msgs/default_topics.h>
#include <mav_msgs/eigen_mav_msgs.h>

#include <ros/ros.h>
#include <std_srvs/Empty.h>
#include <sensor_msgs/Imu.h>
#include <std_msgs/Bool.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <trajectory_msgs/MultiDOFJointTrajectory.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Odometry.h>
//#include <state_graph_builder/graph.h>/
//#include <state_graph_builder/posegraph.h>

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <deque>
#include <string>
#include <algorithm>

struct pose {
    double* header;
    double* x;
    double* y;
    double* theta;
};

struct control_cmd  {
    double* v;
    double* w;
}

class barrier   {
public:
    barrier(int nrovers);
    ~barrier();
private:
    
    // Experiment Parameters
    int nrovers,F,nbad,ngood;
    int[2] badagents;
    int* indices_tot,indices_good,indices_bad;
    
    // Path parameters
    double radius,pi;
    double* circletheta;
    double** tauvector;
    double** y_vector;
    double** prior_y_vector;
    
    // Agent parameters
    double ds,Rs,Rc,el;
    double** Aprox;
    
    // Time parameters
    double t0;
    // Important data structures
    pose state_data,prior_data;
    double** u_data;
    control_cmd input_data;
    
    // ROS node handle
    ros::NodeHandle nh;
    ros::NodeHandle nh_private_;
    
    //msgs
    geometry_msgs::Twist cmd_vel;
    //nav_msgs::Odometry state;
    geometry_msgs::TransformStamped state;
    
    //publisher
    int b;
    string b_str;
    std::string pub_topic;
    ros::Publisher pub;
    ros::Timer pub_timer;
    ros::Timer dis_timer;  // display timer, for debug
    ros::Timer update_param_timer;
    
    // subscriber
    int c;
    string c_str;
    std::string sub_topic;
    ros::Subscriber states_sub;
    
    // Callback functions
    void transformstamped_subCallback(const geometry_msgs::TransformStamped::ConstPtr& msgs);
    void barrierpublisher(const ros::TimerEvent& event);
    void dispublisher(const ros::TimerEvent& event);
    void publishall();
    
    //Mathematical functions
    void update_sparse();
    double Psi_collision_ij(double xi,double yi,double xj,double yj);
    void Psi_gradient_collision(double*& outvector,double xi,double yi,double xj,double yj);
    void velocity_gradient_filtering(vector in_neighbours,vector filtered,vector unfiltered,int index_good);
    void norm_based_filtering(vector<int>& in_neighbours,vector<int>& unfiltered_neighbours,int index_good);
    double Psi_eij(double y_ix,double y_iy,double y_jx,double y_jy,double tauij[2]);
    void Psi_gradient_eij(double*& outvector,double y_ix,double y_iy,double y_jx,double y_jy,double tauij[2]);
    void unicycle_dynamics();
    void calculate_u();
}
#endif /* barrier_h */
