//
//  main.cpp
//  leaderfollower
//
//  Created by Harvard Virgil Humphrey on 2019/7/15.
//  Copyright © 2019 Harvard Virgil Humphrey. All rights reserved.
//

#include "barrier.h"
#include <random>
#include "arrayfunctions.h"
#include "matrixfunctions.h"
#include <vector>
#include "harveymaths.h"

barrier::barrier(int nrovers)
:nh_private_("~")
{
    // Seed random number generator
    srand(time(NULL));
    // Set Parameters
    // Experimental params
    nh_private_.param<int>("nrovers",barrier::nrovers,nrovers);
    nh_private_.param<int>("F",barrier::F,2);
    nh_private_.param<int>("nbad",barrier::nbad,2);
    nh_private_.param<int>("ngood",barrier::ngood,2);
    barrier::indices_tot = new int[nrovers];
    barrier::indices_bad = new int[nbad];
    barrier::indices_good = new int[nrovers-nbad];
    int bob = 0;
    
    // Get a pair of random agent indices between 1 and 6 for malicious agents
    while(1)    {
        if(bob==0)  {
            barrier::indices_bad[0] = (rand()%nrovers)+1;
        }  else {
            barrier""indices_bad[1] = (rand()%nrovers)+1;
        }
        bob++;
        if(barrier::indices_bad[0]!=barrier::indices_bad[1] && bob>=1;)    {
            break;
        }
    }
    int fred = 0;
    // List all agents, and separately, list all non-malicious agents
    for(int i=0;i<nrovers;i++)    {
        barrier::indices_tot[i] = i+1;
        if(barrier::indices_tot[i]!=barrier::indices_bad[0] && barrier::indices_tot[i]!=barrier::indices_bad[1])    {
            barrier::indices_good[fred] = barrier::indices_tot[i];
            fred++;
        }
    }

    
        
    // Path params
    nh_private_.param<double>("radius",radius,1.5);
    nh_private_.param<double>("pi",pi,3.14159265359);
    barrier::circletheta = new double[nrovers];
    zerosdouble(barrier::tauvector,nrovers,2);
    for(int i=0;i<barrier::nrovers;i++) {
        barrier::circletheta[i] = i*(2*pi/barrier::nrovers);
        barrier::tauvector[i][0] = cos(barrier::circletheta[i]);
        barrier::tauvector[i][1] = sin(barrier::circletheta[i]);
    }
    zerosdouble(barrier::y_vector,nrovers,2);
    zerosdouble(barrier::prior_y_vector,nrovers,2);
    
    //Agent params
    nh_private_.param<double>("ds",ds,0.7);
    nh_private_.param<double>("Rs",Rs,0.9);
    nh_private_.param<double>("Rc",Rc,1.1);
    nh_private_.param<double>("el",el,1.0);
    zerosdouble(Aprox,barrier::nrovers,barrier::nrovers);
    
    // Time parameters
    nh_private_.param<double>("t0", barrier::t0, ros::Time().toSec());
    
    // Set pubs and subs
    // Publisher
    barrier::b = 2;
    barrier::b_str = std::to_string(barrier::b);
    barrier::pub_topic = "/R"+barrier::b_str+"/cmd_vel";
    barrier::pub = nh.advertise<geometry_msgs::Twist>(pub_topic, 10);
    barrier::pub_timer = nh.createTimer(ros::Duration(0.02), &barrier::barrierpublisher, this);
    barrier::dis_timer = nh.createTimer(ros::Duration(1), &IO_control_collision::disCallback, this);
    
    // Subscriber
    barrier::c = 2;
    barrier::c_str = std::to_string(barrier::c);
    barrier::sub_topic = "/vicon/R"+barrier::c_str+"/R"+barrier::c_str;
    barrier::states_sub = nh.subscribe(sub_topic,10,&barrier::transformstamped_subCallback,this);
    
    // Initialise structs
    barrier::state_data.x = new double[nrovers];
    barrier::state_data.y = new double[nrovers];
    barrier::state_data.theta = new double[nrovers];
    barrier::prior_data.x = new double[nrovers];
    barrier::prior_data.y = new double[nrovers];
    barrier::prior_data.theta = new double[nrovers];
    barrier::input_data.v = new double[nrovers];
    barrier::input_data.w = new double[nrovers];
    zerosdouble(barrier::u_data,barrier::nrovers,2);
    
    // Populate state structs
    
    
}

void barrier::publishall()  {
    
    // Changes the publisher topic namespace in a loop and publishes multiple messages in rapid succession. Enables quick messaging of all rovers in one line of code
    for(int i=0;i<barrier::nrovers;i++) {
        barrier::b = i+2;
        barrier::b_str = std::to_string(barrier::b);
        barrier::pub_topic ="/R"+barrier::b_str+"/cmd_vel";
        barrier::pub = nhadvertise<geometry_msgs::Twist>(pub_topic,10);
        barrier::cmd_vel.linear.x = barrier::input_data.v[i];
        barrier::cmd_vel.angular.z = barrier::input_data.w[i];
        pub.publish(cmd_vel);
    }
}

barrier::~barrier() {
    barrier::cmd_vel.linear.x = 0;
    barrier::cmd_vel.angular.z = 0;
    for(int i=0;i<barrier::nrovers;i++) {
        barrier::b = i+2;
        barrier::b_str = std::to_string(barrier::b);
        barrier::pub_topic ="/R"+barrier::b_str+"/cmd_vel";
        barrier::pub = nhadvertise<geometry_msgs::Twist>(pub_topic,10);
        pub.publish(barrier::cmd_vel);
    }
    ros::shutdown();
    
    delete[] barrier::state_data.x;
    delete[] barrier::state_data.y;
    delete[] barrier::state_data.theta;
    delete[] barrier::prior_data.x;
    delete[] barrier::prior_data.y;
    delete[] barrier::prior_data.theta;
    delete[] barrier::input_data.v;
    delete[] barrier::input_data.w;
    delete[] barrier::circletheta;
    deletematrixdouble(barrier::tauvector,barrier::nrovers,2);
    deletematrixdouble(barrier::y_vector,barrier::nrovers,2)l
    deletematrixdouble(barrier::prior_y_vector,barrier::nrovers,2);
    deletematrixdouble(barrier::Aprox,barrier::nrovers,barrier::nrovers);
    deletematrixdouble(barrier::u_data,barrier::nrovers,2);
}

void barrier::transformstamped_subCallback(const geometry_msgs::TransformStamped::ConstPtr& msgs)   {
    for(int i=0;i<barrier::nrovers;i++) {
        double head1,head2,head3,head4;
        // Time reading
        barrier::state_data.header[i] = msgs->header;
        barrier::state_data.x[i] = msgs->transform.translation.x;
        barrier::state_data.y[i] = msgs->transform.translation.y;
        
        // Orientation calculation
        barrier::state_data.theta[i] = msgs->transform.translation.z;
        head1 = msgs->transform.rotation.x;
        head2 = msgs->transform.rotation.y;
        head3 = msgs->transform.rotation.z;
        head4 = msgs->transform.rotation.w;
        tf::Quaternion q(head1,head2,head3,head4);
        tf::Matrix3x3 m(q);
        double roll,pitch,yaw;
        m.getRPY(roll,pitch,yaw);
        barrier::theta[i] = fmod(yaw+2*M_PI,2*M_PI);
    }
}
void barrier::update_sparse()   {
    // Identifying all neighbouring agents with respect to agent i and the safety distsance
    double stateix,stateiy,statejx,statejy,diffx,diffy,dist
    for(int i=0;i<barrier::nrovers;i++) {
        stateix = barrier::state_data.x[i];
        stateiy = barrier::state_data.y[i];
        for(int j=0;j<barrier::nrovers;j++) {
            statejx = barrier::state_data.x[j];
            statejy = barrier::state_data.y[j];
            diffx = abs(stateix-statejx);
            diffy = abs(stateiy-statejy);
            dist = sqrt(diffx*diffx+diffy*diffy);
            if(dist<=barrier::el)   {
                barrier::Aprox[i][j] = 1;
            }
        }
    }
}

double barrier::Psi_collision_ij(double xi,double yi,double xj,double yj)    {
    double diffX = xi-xj;
    double diffY = yi-yj;
    double norm = sqrt(diffX*diffX+diffY*diffY);
    double mu2 = 10000.0;
    double outscalar;
    if(norm<=barrier::Rs)    {
        outscalar = (norm-barrier::Rs)*(norm-barrier::Rs)/(norm-barrier::Rs+(barrier::ds-barrier::Rs)*(barrier::ds-barrier::Rs)/mu2);
    }
    else if(norm<barrier::ds)  {
        outscalar = mu2;
    }   else {
        outscalar = 0;
    }
    return outscalar;
}

void barrier::Psi_gradient_collision(double*& outvector,double xi,double yi,double xj,double yj)   {
    double h = 0.001;
    outvector[0] = (Psi_collision_ij(xi+h,yi,xj,yj)-Psi_collision_ij(xi-h,yi,xj,yj))/(2*h);
    outvector[1] = (Psi_collision_ij(xi,yi+h,xj,yj)-Psi_collision_ij(xi,yi-h,xj,yj))/(2*h);
}

double barrier::Psi_eij(double y_ix,double y_iy,double y_jx,double y_jy,double tauij[2])   {
    double mu1=1000;
    double yij[2];
    yij[0] = y_ix-y_jx;
    yij[1] = y_iy-y_jy;
    
    double rshat = barrier::el-sqrt(tauij[0]*tauij[0]+tauij[1]*tauij[1]);
    double norm_y = sqrt(yij[0]*yij[0]+yij[1]*yij[1]);
    double outscalar = (norm_y*norm_y)/(rshat-norm_y+(rshat*rshat)/mu1);
    return outscalar;
}

void barrier::Psi_gradient_eij(double*& outvector,double y_ix,double y_iy,double y_jx,double y_jy,double tauij[2])  {
    double h = 0.001;
    outvector[0] = (barrier::Psi_eij(y_ix+h,y_iy,y_jx,y_jy,tauij)-Psi_eij(y_ix-h,y_iy,y_jx,y_jy,tauij))/(2*h);
    outvector[1] = (barrier::Psi_eij(y_ix,y_iy+h,y_jx,y_jy,tauij)-Psi_eij(y_ix,y_iy-h,y_jx,y_jy,tauij))/(2*h);
}

void barrier::velocity_gradient_filtering(vector in_neighbours,vector filtered_neighbours,vector unfiltered_neighbours,int index_good) {
    bool empty = in_neighbours.empty();
    double** yvel_times_T;
    double** gradient_vector;
    double** yij;
    double** velocity_times_gradient;
    double** filtered_neighbours;
    double** unfiltered_neighbours;
    int size = in_neighbours.size();
    zerosdouble(yij,barrier::ngood,2);
    zerosdouble(yvel_times_T,barrier::nrovers,2);
    zerosdouble(gradient_vector,barrier::ngood,2);
    zerosdouble(velocity_times_gradient,barrier::ngood,2);
    if(empty==false) {
        // Technically yvel would be calculated as yvel = 1/T(state_vector - prior_state_vector). But since we're just sorting by relative magnitude, and 1/T is a positive constant, the 1/T term can be omitted without changing the order of the values w.r.t the relative magnitude between vectors.
        for(int i=0;i<barrier::nrovers;i++) {
            for(int j=0;j<2;j++)    {
                yvel_times_T[i][j] = barrier::y_vector[i][j]-barrier::prior_y_vector[i][j];
            }
        }
        double yidot[2];
        yidot[0] = yvel_times_T[index_good][0];
        yidot[1] = yvel_times_T[index_good][1];
        int jj;
        double tauij[2];
        for(int i=0;i<size;i++)   {
            jj = in_neighbours.at(i);
            jj = jj-1;
            tauij[0] = tauvector[index_good][0]-tauvector[jj][0];
            tauij[1] = tauvector[index_good][1]-tauvector[jj][1];
            Psi_gradient_eij(gradient_vector[i],barrier::y_vector[index_good][0],barrier::y_vector[index_good][1],barrier::y_vector[jj][0],barrier::y_vector[jj][1],tauij);
            yij[i][0] = yidot[0]-yvel_times_T[jj][0];
            yij[i][1] = yidot[1]-yvel_times_T[jj][1];
        }
        
        // Create the list
        for(int i=0;i<size;i++)   {
            velocity_times_gradient[i][0] = yij[i][0]*gradient_vector[i][0]+yij[i][1]*gradient_vector[i][1];
            velocity_times_gradient[i][1] = barrier::indices_good[i];
        }
        
        // Sort the list
        sortrowsdescend(velocity_times_gradient,0,size,2);
        
        // Choose which list to go with
        if(barrier::F<=size)  {
            
            // Need to empty the vectors here first
            for(int i=0;i<barrier::F;i++)   {
                unfiltered_neighbours.push_back(velocity_times[i][1]);
            }
            for(int i=F;i<size;i++)   {
                filtered_neighbours.push_back(velocity_times[i][1]);
            }
        }   else    {
            
            // Need to empty vectors here first
            for(int i=0;i<size;i++) {
                filtered_neighbours.push_back(velocity_times[i][1]);
            }
        }
    }   else    {
        // Empty the vectors
    }
    deletematrixdouble(yij,size,2);
    deletematrixdouble(yvel_times_T,barrier::nrovers,2);
    deletematrixdouble(gradient_vector,barrier::size,2);
    deletematrixdouble(velocity_times_gradient,barrier::size,2);
}

void barrier::norm_based_filtering(vector<int>& in_neighbours,vector<int>& unfiltered_neighbours,int index_good) {
    double y_ix = barrier::y_vector[index_good-1][0];
    double y_iy = barrier::y_vector[index_good-1][0];
    double diffX,diffY,y_jx,y_jy;
    double** norm_diff_vector;
    int size = in_neighbours.size();
    zerosdouble(norm_diff_vector,size,2);
    int jj;
    for(int i=0;i<size;i++) {
        jj = in_neighbours.at(i);
        jj = jj-1;
        y_jx = barrier::y_vector[jj][0];
        y_jy = barrier::y_vector[jj][1];
        diffX = y_ix-y_jx;
        diffY = y_iy-y_jy;
        norm_diff_vector[i][0] = sqrt(diffX*diffX+diffY*diffY);
        norm_diff_vector[i][1] = jj+1;
    }
    sortrowsdescend(norm_diff_vector,1,size,2);
    
    for(int i=0;i<size-barrier::F;i++)  {
        unfiltered_neighbours.push_back(norm_diff_vector[i+F-1][1]);
    }
    deletematrixdouble(norm_diff_vector,size,2);
}

void barrier::calculate_u() {
    int index;
    double* outvector = new double[2];
    double* Psi_collision_sum = new double[2];
    zerosarraydouble(Psi_collision_sum,2);
    // Get the right vectors
    for(int i=0;i<barrier::nrovers;i++) {
        barrier::y_vector[i][0] = barrier::state_data.x[i]-barrier::tauvector[i][0];
        barrier::y_vector[i][1] = barrier::state_data.y[i]-barrier::tauvector[i][1];
        barrier::prior_y_vector[i][0] = barrier::prior_data.x[i]-barrier::tauvector[i][0];
        barrier::prior_y_vector[i][1] = barrier::prior_data.y[i]-barrier::tauvector[i][1];
    }
    
    // We address the malicious agents first
    for(int i=0;i<barrier::nbad;i++)    {
        index = barrier::indices_bad[i];
        barrier::input_data.v[index-1] = 0;
        barrier::input_data.w[index-1] = 0;
    }
    
    // Now we do the good agents
    std::vector<int> in_neighbours;
    std::vector<int> filtered_neighbours;
    std::vector<int> unfiltered_neighbours;
    int index_good;
    int size;
    int indexj;
    double statex,statey,statejx,statejy;
    for(int i=0;i<barrier::ngood;i++)   {
        insertzerosdouble(Psi_collision_sum,2);
        index_good = barrier::indices_good[i];
        index_good--;
        statex = barrier::state_data.x[index_good];
        statey = barrier::state_data.y[index_good];
        
        // Find the elements of Aprox that are 1 for each agent row
        for(int j=0;j<barrier::nrovers;j++) {
            if(barrier::Aprox[index_good][j]==1)  {
                in_neighbours.push_back(j+1);
            }
        }
        
        // Calculate the collision avoidance barrier function. The collision avoidance barrier is unfiltered--it avoids all agents within the correct proximity.
        size = in_neighbours.size();
        if(in_neighbours.empty()==0)    {
            for(int j=0;j<size;j++) {
                indexj = in_neighbours.at(j);
                indexj = indexj-1;
                statejx = barrier::state_data.x[indexj];
                statejy = barrier::state_data.y[indexj];
                Psi_gradient_collision(outvector,statex,statey,statejx,statejy);
                Psi_collision_sum[0]+=outvector[0];
                Psi_collision_sum[1]+=outvector[1];
            }
        }
        
        // Create the filtered in-neighbour list
        //barrier::velocity_gradient_filtering(in_neighbours,filtered_neighbours,unfiltered_neighbours,i);
        barrier::norm_based_filtering(in_neighbours,unfiltered_neighbours,i);
        double yi_x,yi_y,yj_x,yj_y;
        yi_x = barrier::y_vector[i][0];
        yi_y = barrier::y_vector[i][1];
        int index;
        double tauij[2];
        double* Psi_gradient_sum = new double[2];
        zerosarraydouble(Psi_gradient_sum,2);
        if(unfiltered_neighbours.empty()==false)    {
            for(int jj=0;jj<unfiltered_neighbours.size();jj++)   {
                index = unfiltered_neighbours.at(jj);
                index = index-1;
                yj_x = barrier::y_vector[index][0];
                yj_y = barrier::y_vector[index][0];
                tauij[0] = barrier::tauvector[i][0]-barrier::tauvector[index][0];
                tauij[1] = barrier::tauvector[i][1]-barrier::tauvector[index][1];
                Psi_gradient_eij(outvector,yi_x,yi_y,yj_x,yj_y,tauij);
                Psi_gradient_sum[0]+=outvector[0];
                Psi_gradient_sum[1]+=outvector[1];
            }
        }
        double gain = 1.0;
        barrier::u_data[i][0] = -(Psi_gradient_sum[0]+Psi_collision_sum[0])*gain;
        barrier::u_data[i][1] = -(Psi_gradient_sum[1]+Psi_collision_sum[1])*gain;\
        delete[] Psi_gradient_sum;
    }
    delete[] outvector;
    delete[] Psi_collision_sum;
}

void barrier::unicycle_dynamics()   {
    
    double umin_linear;
    double umin_angular;
    double umax_linear = 0.5;
    double umax_angular = 3.0;
    double linear_error,linear_temp,turn_error,turn_temp,turn_req,heading_req;
    double gain_linear = 0.1525;
    double gain_angular = 0.3050;
    for(int i=0;i<barrier::nrovers;i++)    {
        heading_req = atan2(barrier::u_data[i][1],barrier::u_data[i][0]);
        turn_req = heading_req-barrier::state_data.theta[i];
        turn_error = wrap2pi(turn_req);
        linear_error = sqrt(barrier::u_data[i][0]*barrier::u_data[i][0]+barrier::u_data[i][1]*barrier::u_data[i][1]);
        
        linear_temp = umin_linear+linear_error*gain_linear;
        if(linear_temp>umax_linear) {
            linear_temp = umax_linear;
        }
        if(turn_error>0.05) {
            linear_temp = 0.0;
        }
        if(linear_temp=0.0) {
            umin_angular = 1.0;
        }   else    {
            umin_angular = 0.6;
        }
        
        if(turn_error>0)    {
            turn_temp = gain_angular*turn_error+umin_angular;
            if(turn_temp>umax_angular)  {
                turn_temp = umax_angular;
            }
        }
        else if(turn_error<0)   {
            turn_temp = gain_angular*turn_error-umin_angular;
            if(turn_temp<-umax_angular) {
                turn_temp = -umax_angular;
            }
        }   else    {
            turn_temp = 0;
        }
        
        
        // Assign to the publishing message
        barrier::input_data.v[i] = linear_temp;
        barrier::input_data.w[i] = turn_temp;
    }
    
}

void barrier::barrierpublisher(const ros::TimeEvent& event) {
    barrier::calculate_u();
    barrier::unicycle_dynamics();
    barrier::publishall();
}

int main(int argc, const char * argv[]) {
    ros::init(argc,argv,"Resilient_Barrier");
    barrier Resilient_Barrier(6);
    ros::spin();
    return 0;
}
