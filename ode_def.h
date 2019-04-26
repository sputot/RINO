/* ============================================================================
 File   : ode_def.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 The place to declare the systems of ODEs or DDEs on which he/she wants to perform reachability
 Setup of dimension of the system (sysdim), initial conditions and parameters is in ode_def.cpp
 ============================================================================ */
#ifndef ODE_DEF_H
#define ODE_DEF_H

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"

using namespace std;

#define innerapprox 1

extern int systype;    // systype = 0 (ODE) or 1 (DDE) -- initialized in main.cpp / command line
extern int syschoice;  // to choose among the predefined systems of ODE or DDE -- initialized in main.cpp / command line

// dimension of the system of ODE/DDE to analyze:
extern int sysdim;
// dimension of uncertain input
extern int jacdim; // Jacobian will be dimension sysdim * jacdim
extern int sysdim_params;  // dimension of the vector of parameters params - Boost
extern vector<AAF> params;      // params of the ODE (nondeterministic disturbances)
extern vector<AAF> inputs;   // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
extern vector<AAF> center_inputs;
extern vector<interval> eps;

// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
extern vector<double> t_print; // times where results are stored
extern int current_subdiv;
extern int current_iteration;

// for robust inner-approximations
extern int uncontrolled;  // number of uncontrolled parameters (forall params)
extern int controlled;  // number of controlled parameters (forall params)
extern vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
extern vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust inner-approx)
extern int variable;  // number of non constant parameters
extern vector<bool> is_variable; // for each parameter, constant or variable


void define_system_dim();  // define the dimensions of your system (ODE or DDE)




// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for DDEs : functions that initialize the DDE on first time period
vector<T<AAF>> Initfunc(const T<AAF>& t, vector<AAF> &beta);
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta);

// defining analytical solution if any for comparison
vector <interval> AnalyticalSol(double t, vector<AAF> &beta, double d0);


// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide);



// define here  your ODE system yp = \dot y = f(y)
class OdeFunc {
public:
    template <class C>
      void operator()(vector<C> &yp, vector<C> y) {
          
          
          if (syschoice == 1) // running example
             yp[0] = (1+params[0]*y[0])*(1+params[0]*y[0]);
          else if (syschoice == 2)
          // Bruss corresponds to the Brusselator system  
          {
              yp[0] = 1 - (params[1]+1)*y[0] + params[0]*y[0]*y[0]*y[1];
              yp[1] = params[1]*y[0] - params[0]*y[0]*y[0]*y[1];
          }
          else if (syschoice == 3)
          /******************* ballistic ****************************/
          {
           double g = 9.81; // gravity in m/s^2
           double rho = 1204.4; // air density in g/m3
           double a = 0.000126677; // cross section du bullet d=12.7mm, cross section=pi R^2
           double d = 0.45; // drag coefficient
           double m = 14.3; // mass du bullet in g
           yp[0] = - g*sin(y[1])-rho*y[0]*y[0]*a*d/(2.0*m); // velocity v
           yp[1] = - g*cos(y[1])/y[0]; // angle gamma with respect to the x axis
           yp[2] = y[0]*cos(y[1]); // position x
           yp[3] = y[0]*sin(y[1]); // position y
          }
          else if (syschoice == 4)
          /**** ballistic linearise ... *****/
          {
           double g = 9.81; // gravity in m/s^2
           double rho = 1204.4; // air density in g/m3
           double a = 0.000126677; // cross section du bullet d=12.7mm, cross section=pi R^2
           double d = 0.45; // drag coefficient
           double m = 14.3; // mass du bullet in g
         //  yp[0] = - g*(y[1])-rho*y[0]*y[0]*a*d/(2.0*m); // velocity v
              yp[0] = - g*(y[1])-rho*y[0]*y[0]*a*d/(2.0*y[4]); // velocity v
           yp[1] = - g*(1-y[1]*y[1]/2)/y[0]; // angle gamma with respect to the x axis
           yp[2] = y[0]*(1-y[1]*y[1]/2); // position x
           yp[3] = y[0]*(y[1]); // position y
              yp[4] = 0;
          }
          /******************* end ballistic ****************************/
          else if (syschoice == 5)  // self-driving car
          {
              yp[0] = y[1];
              yp[1] = -params[0] *(y[0] - 1.0) - params[1]*y[1];   // pr = 1 is the reference position
          }
        /*  else if (syschoice == 6)  // self-driving car (version obsolete)
          {
              yp[0] = y[1];
              yp[1] = -inputs[2] *(y[0] - 1.0) - inputs[3]*y[1];   // pr = 1 is the reference position
          }*/
          else if (syschoice == 6)  // self-driving car with constant parameters
          {
              yp[0] = y[1];
              yp[1] = -y[2] *(y[0] - 1.0) - y[3]*y[1];   // pr = 1 is the reference position
              yp[2] = 0; // constant parameter Kp
              yp[3] = 0; // constant parameter Kd
          }
          else if (syschoice == 7)  // self-driving car with time-varying parameters
          {
              yp[0] = y[1];
              yp[1] = -y[2] *(y[0] - 1.0) - y[3]*y[1];   // pr = 1 is the reference position
              yp[2] = interval(-2.,2.); //AAF(interval(-2.,2.)) + y[4];  // parameter Kp
              yp[3] = interval(-2.,2.); // 0;  // parameter Kd
          }
          else if (syschoice == 8) // -f for viability problem
          {
              yp[0] = y[1] - y[2];
              yp[1] = - y[0] + y[1]*y[1]*y[1] - y[3];
              yp[2] = 0;  // uncontrollable parameter
              yp[3] = 0; // uncontrollable parameter
           /*   yp[4] = y[5] - y[6];
              yp[5] = - y[4] + y[5]*y[5]*y[5] - y[7];
              yp[6] = 0;  // uncontrollable parameter
              yp[7] = 0; // uncontrollable parameter
              yp[8] = y[9] - y[10];
              yp[9] = - y[8] + y[9]*y[9]*y[9] - y[11];
              yp[10] = 0;  // uncontrollable parameter
              yp[11] = 0; // uncontrollable parameter
              yp[12] = y[13] - y[14];
              yp[13] = - y[12] + y[13]*y[13]*y[13] - y[15];
              yp[14] = 0;  // uncontrollable parameter
              yp[15] = 0; // uncontrollable parameter */
          }
          else if (syschoice == 9) // Crazyflie attitude control problem : version gitlab F. Djeumou
          {
              static const double thrust = 35000;
              
              static const double roll_sp =10.0;
              static const double pitch_sp = 0.0;
              static const double yaw_sp = 0.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const AAF Ct = interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_r = 150.0;
              static const double Kp_p = 150.0;
              static const double Kp_y = 100.0;
              
              static const double Kd_r = 120;
              static const double Kd_p = 120;
              static const double Kd_y = 120;
              
              static const double Ki_r = 30.0;
              static const double Ki_p = 30.0;
              static const double Ki_y = 10.0;
              
              // static const double Kp_r = 300.0;
              // static const double Kp_p = 300.0;
              // static const double Kp_y = 150.0;
              
              // static const double Kd_r = 100;
              // static const double Kd_p = 100;
              // static const double Kd_y = 100;
              
              // static const double Ki_r = 60.0;
              // static const double Ki_p = 60.0;
              // static const double Ki_y = 30.0;
              
              // static const AAF Ip_p = interval(-9.884e-4 , 0.0);
              // static const AAF Ip_pq = interval(-4.884e-2 , 0.0);
              // static const AAF Ip_pr = interval(0.0 , 1.643e-2);
              // static const AAF Ip_q = interval(-1.104e-1 , 0.0);
              // static const AAF Ip_qr = interval(-1.048 , -0.7571);
              // static const AAF Ip_r = interval(0.0 , 0.1115);
              
              // static const AAF Iq_p = interval(0.0 , 4.673e-4);
              // static const AAF Iq_pq = interval(0.0 , 0.112);
              // static const AAF Iq_pr = interval(0.7681 , 1.027);
              // static const AAF Iq_q = interval(0.0 , 2.16e-3);
              // static const AAF Iq_qr = interval(0.0 , 4.455e-2);
              // static const AAF Iq_r = interval(-4.887e-2 , 0.0);
              
              // static const AAF Ir_p = interval(-3.133e-2 , 0.0);
              // static const AAF Ir_pq = interval(-1.629e-2 , 2.08e-3);
              // static const AAF Ir_pr = interval(-0.11 , 0.0);
              // static const AAF Ir_q = interval(0.0 , 3.107e-2);
              // static const AAF Ir_qr = interval(0.0 , 4.455e-2);
              // static const AAF Ir_r = interval(0.0 , 2.797e-4);
              
              // static const AAF Im_xx = interval(60398.23 , 71885.73);
              // static const AAF Im_xy = interval(-2886.22 , 0.0);
              // static const AAF Im_xz = interval(-1313.91 , 0.0);
              // static const AAF Im_yy = interval(60446.53 , 70125.47);
              // static const AAF Im_yz = interval(-3667.42 , 0.0);
              // static const AAF Im_zz = interval(34337.83 , 34712.03);
              
              static const AAF Ip_p = 0.0;
              static const AAF Ip_pq = 0.0;
              static const AAF Ip_pr = 0.0;
              static const AAF Ip_q = 0.0;
              static const AAF Ip_qr = interval(-1.04880447793, -1.03580464787);
              static const AAF Ip_r = 0.0;
              
              static const AAF Iq_p = 0.0;
              static const AAF Iq_pq = 0.0;
              static const AAF Iq_pr = interval(1.03470095927, 1.04749270535);
              static const AAF Iq_q = 0.0;
              static const AAF Iq_qr = 0.0;
              static const AAF Iq_r = 0.0;
              
              static const AAF Ir_p = 0.0;
              static const AAF Ir_pq = interval(-0.0162919189567, -0.0120891632629);
              static const AAF Ir_pr = 0.0;
              static const AAF Ir_q = 0.0;
              static const AAF Ir_qr = 0.0;
              static const AAF Ir_r = 0.0;
              
              static const AAF Im_xx = interval(71484.0524534, 71885.7226787);
              static const AAF Im_xy = 0.0;
              static const AAF Im_xz = 0.0;
              static const AAF Im_yy = interval(69441.6509547, 69834.7034512);
              static const AAF Im_yz = 0.0;
              static const AAF Im_zz = interval(34492.4780616, 34712.0265858);
              
              auto cosy0 = cos(y[0]*M_PI/180.0);
              auto siny0 = sin(y[0]*M_PI/180.0);
              auto tany1 = tan(y[1]*M_PI/180.0);
              // Rotationnal kinematic euler angle version
              yp[0] = y[3] + (y[5]*cosy0 + y[4]*siny0)*tany1;
              yp[1] = y[4]*cosy0 - y[5]*siny0;
              yp[2] = (y[5]*cosy0 + y[4]*siny0)/cos(y[1]*M_PI/180.0);
              
              yp[6] = roll_sp - y[0];
              yp[7] = pitch_sp - y[1];
              yp[8] = yaw_sp - y[2];
              
              auto cmd_r = y[6]*Ki_r + yp[6]*Kp_r - yp[0]*Kd_r;
              auto cmd_p = y[7]*Ki_p + yp[7]*Kp_p - yp[1]*Kd_p;
              auto cmd_y = y[8]*Ki_y + yp[8]*Kp_y - yp[2]*Kd_y;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              
              yp[3] = Ip_p * y[3] * y[3] + Ip_pq * y[3] * y[4] + Ip_pr * y[3] * y[5] + Ip_q * y[4] * y[4] + Ip_qr * y[4] * y[5] + Ip_r * y[5] * y[5] + Im_xx * Mx + Im_xy * My + Im_xz * Mz;
              yp[4] = Iq_p * y[3] * y[3] + Iq_pq * y[3] * y[4] + Iq_pr * y[3] * y[5] + Iq_q * y[4] * y[4] + Iq_qr * y[4] * y[5] + Iq_r * y[5] * y[5] + Im_xy * Mx + Im_yy * My + Im_yz * Mz;
              yp[5] = Ir_p * y[3] * y[3] + Ir_pq * y[3] * y[4] + Ir_pr * y[3] * y[5] + Ir_q * y[4] * y[4] + Ir_qr * y[4] * y[5] + Ir_r * y[5] * y[5] + Im_xz * Mx + Im_yz * My + Im_zz * Mz;
              
          }
          else if (syschoice == 10) // crazyflie inner cascade loop --> No uncertaincy are used here
          {
              // The formula here are obtained using symbolic calculus in matlab
              // And by solving a DAE (acttually it's linear here) to an ODE
              
              // Should be in radian
              static const double roll_sp = 10.0 * M_PI/180.0;
              static const double pitch_sp = 0.0;
              static const double yaw_sp = 0.0;//7.0 * M_PI/180.0;
              static const double thrust = 50000;
              
              auto cosy0 = cos(y[0]);
              auto siny0 = sin(y[0]);
              auto tany1 = tan(y[1]);
              
              yp[0] = y[3] + (y[5]*cosy0 + y[4]*siny0)*tany1;
              yp[1] = y[4]*cosy0 - y[5]*siny0;
              yp[2] = (y[5]*cosy0 + y[4]*siny0)/cos(y[1]);
              
              
              // integral of p , q and r
              yp[6] = y[3];
              yp[7] = y[4];
              yp[8] = y[5];
              
              // integral of error in roll , pitch and yaw
              yp[9] = y[0];
              yp[10] = y[1];
              yp[11] = y[2];
              
              // Double integrale of error in roll , pitch and yaw
              yp[12] = y[9];
              yp[13] = y[10];
              yp[14] = y[11];
              
              // Time
              yp[15] = 1.0;
              
              auto sp_yaw_d = 0.0;
              auto sp_yaw_i = yaw_sp * y[15];
              auto sp_yaw_ii = sp_yaw_i * y[15] * 0.5;
              
              auto sp_roll_d = 0.0;
              auto sp_roll_i = roll_sp * y[15];
              auto sp_roll_ii = sp_roll_i * y[15] * 0.5;
              
              auto sp_pitch_d = 0.0;
              auto sp_pitch_i = pitch_sp * y[15];
              auto sp_pitch_ii = sp_pitch_i * y[15] * 0.5;
              
              // auto k1_p = -0.02906;
              auto k1_p = - 4.191e-7*thrust - 0.003916;
              auto k2_p = 0.0003042*yaw_sp - 5.03e-5*y[5] - 0.0003042*y[2] + 7.0e-6*sp_yaw_ii - 7.0e-6*y[14] + 1.76e-5*sp_yaw_d - 7.0e-6*y[8] + 9.229e-5*sp_yaw_i - 1.76e-5*yp[2] - 9.229e-5*y[11];
              auto k1_q = 0.0003027*yaw_sp - 5.004e-5*y[5] - 0.0003027*y[2] + 6.964e-6*sp_yaw_ii - 6.964e-6*y[14] + 1.752e-5*sp_yaw_d - 6.964e-6*y[8] + 9.183e-5*sp_yaw_i - 1.752e-5*yp[2] - 9.183e-5*y[11];
              auto k2_q = - 4.17e-7*thrust - 0.003896;
              
              auto denom  = (1.0 - k1_p) * (1.0 - k2_q) - k1_q*k2_p;
              //auto k3_p = 17.53*roll_sp - 17.53*y[0] - 2.906*y[3] - 0.0003018*y[5]*yp[1] + 0.01062*y[1]*sp_yaw_d + 0.001825*yaw_sp*yp[1] + 0.00176*y[4]*sp_yaw_d + 0.0003018*y[5]*sp_pitch_d - 0.001825*yaw_sp*sp_pitch_d - 0.01062*pitch_sp*sp_yaw_d - 0.004221*y[1]*y[8] - 0.07545*y[5]*y[10] - 0.0007*y[4]*y[8] - 0.01006*y[5]*y[7] + 0.05565*y[1]*sp_yaw_i + 0.4563*yaw_sp*y[10] + 0.009229*y[4]*sp_yaw_i + 0.06085*yaw_sp*y[7] - 0.01062*y[1]*yp[2] - 0.001825*y[2]*yp[1] - 0.00176*y[4]*yp[2] + 0.07545*y[5]*sp_pitch_i + 0.004221*pitch_sp*y[8] - 0.4563*yaw_sp*sp_pitch_i - 0.05565*pitch_sp*sp_yaw_i + 0.01062*pitch_sp*yp[2] + 0.001825*y[2]*sp_pitch_d - 0.05565*y[1]*y[11] - 0.4563*y[2]*y[10] - 0.009229*y[4]*y[11] - 0.06085*y[2]*y[7] + 0.05565*pitch_sp*y[11] + 0.4563*y[2]*sp_pitch_i - 0.03033*y[1]*y[5] - 0.7657*y[4]*y[5] + 0.1835*y[1]*yaw_sp + 0.03042*y[4]*yaw_sp + 0.03033*y[5]*pitch_sp - 0.1835*yaw_sp*pitch_sp - 0.1835*y[1]*y[2] - 0.03042*y[4]*y[2] + 0.1835*pitch_sp*y[2] + 0.0042*y[13]*sp_yaw_ii - 17.44*y[12] + 17.44*sp_roll_ii - 0.0042*sp_yaw_ii*sp_pitch_ii - 0.0042*y[13]*y[14] + 0.0042*sp_pitch_ii*y[14] + 0.01056*y[13]*sp_yaw_d + 4.2e-5*sp_yaw_ii*yp[1] - 0.1744*yp[0] + 0.1744*sp_roll_d - 5.813*y[6] - 4.2e-5*sp_yaw_ii*sp_pitch_d - 0.01056*sp_pitch_ii*sp_yaw_d - 0.0042*y[13]*y[8] + 0.05538*y[13]*sp_yaw_i + 0.0105*sp_yaw_ii*y[10] - 43.6*y[9] + 0.0014*sp_yaw_ii*y[7] - 0.01056*y[13]*yp[2] - 4.2e-5*y[14]*yp[1] + 43.6*sp_roll_i + 0.0042*sp_pitch_ii*y[8] - 0.0105*sp_yaw_ii*sp_pitch_i - 0.05538*sp_pitch_ii*sp_yaw_i + 0.01056*sp_pitch_ii*yp[2] + 4.2e-5*y[14]*sp_pitch_d - 0.05538*y[13]*y[11] - 0.0105*y[14]*y[10] - 0.0014*y[14]*y[7] + 0.05538*sp_pitch_ii*y[11] + 0.0105*y[14]*sp_pitch_i + 0.0001056*yp[1]*sp_yaw_d - 0.0001056*sp_yaw_d*sp_pitch_d - 4.2e-5*yp[1]*y[8] + 0.0005538*yp[1]*sp_yaw_i + 0.02641*sp_yaw_d*y[10] + 0.003521*sp_yaw_d*y[7] - 0.0001056*yp[1]*yp[2] + 4.2e-5*sp_pitch_d*y[8] - 0.02641*sp_yaw_d*sp_pitch_i - 0.0005538*sp_pitch_d*sp_yaw_i - 0.03018*y[5]*y[13] - 0.0105*y[10]*y[8] + 0.0001056*sp_pitch_d*yp[2] - 0.0014*y[7]*y[8] + 0.004221*y[1]*sp_yaw_ii + 0.1825*yaw_sp*y[13] + 0.1384*y[10]*sp_yaw_i + 0.0007*y[4]*sp_yaw_ii + 0.01846*y[7]*sp_yaw_i - 0.0005538*yp[1]*y[11] - 0.02641*y[10]*yp[2] - 0.003521*y[7]*yp[2] + 0.03018*y[5]*sp_pitch_ii + 0.0105*y[8]*sp_pitch_i - 0.1825*yaw_sp*sp_pitch_ii - 0.004221*pitch_sp*sp_yaw_ii - 0.1384*sp_yaw_i*sp_pitch_i + 0.0005538*sp_pitch_d*y[11] + 0.02641*sp_pitch_i*yp[2] - 0.004221*y[1]*y[14] - 0.1825*y[2]*y[13] - 0.1384*y[10]*y[11] - 0.0007*y[4]*y[14] - 0.01846*y[7]*y[11] + 0.004221*pitch_sp*y[14] + 0.1825*y[2]*sp_pitch_ii + 0.1384*sp_pitch_i*y[11];
              auto k3_p = 2.361*roll_sp - 2.361*y[0] - 0.3916*y[3] - 0.0003018*y[5]*yp[1] + 0.01062*y[1]*sp_yaw_d + 0.001825*yaw_sp*yp[1] + 0.00176*y[4]*sp_yaw_d + 0.0003018*y[5]*sp_pitch_d - 2.515e-6*thrust*yp[0] - 0.001825*yaw_sp*sp_pitch_d - 0.01062*pitch_sp*sp_yaw_d - 0.004221*y[1]*y[8] - 0.07545*y[5]*y[10] + 2.515e-6*thrust*sp_roll_d - 0.0007*y[4]*y[8] - 0.01006*y[5]*y[7] + 0.05565*y[1]*sp_yaw_i + 0.4563*yaw_sp*y[10] - 8.383e-5*thrust*y[6] + 0.009229*y[4]*sp_yaw_i + 0.06085*yaw_sp*y[7] - 0.01062*y[1]*yp[2] - 0.001825*y[2]*yp[1] - 0.00176*y[4]*yp[2] + 0.07545*y[5]*sp_pitch_i + 0.004221*pitch_sp*y[8] - 0.0006287*thrust*y[9] - 0.4563*yaw_sp*sp_pitch_i - 0.05565*pitch_sp*sp_yaw_i + 0.01062*pitch_sp*yp[2] + 0.001825*y[2]*sp_pitch_d + 0.0006287*thrust*sp_roll_i - 0.05565*y[1]*y[11] - 0.4563*y[2]*y[10] - 0.009229*y[4]*y[11] - 0.06085*y[2]*y[7] + 0.05565*pitch_sp*y[11] + 0.4563*y[2]*sp_pitch_i - 0.03033*y[1]*y[5] - 0.7657*y[4]*y[5] + 0.1835*y[1]*yaw_sp - 4.191e-5*y[3]*thrust + 0.03042*y[4]*yaw_sp + 0.03033*y[5]*pitch_sp - 0.0002527*y[0]*thrust - 0.1835*yaw_sp*pitch_sp + 0.0002527*roll_sp*thrust - 0.1835*y[1]*y[2] - 0.03042*y[4]*y[2] + 0.1835*pitch_sp*y[2] + 0.0042*y[13]*sp_yaw_ii - 2.349*y[12] + 2.349*sp_roll_ii - 0.0042*sp_yaw_ii*sp_pitch_ii - 0.0042*y[13]*y[14] + 0.0042*sp_pitch_ii*y[14] + 0.01056*y[13]*sp_yaw_d + 4.2e-5*sp_yaw_ii*yp[1] - 0.02349*yp[0] + 0.02349*sp_roll_d - 0.7831*y[6] - 4.2e-5*sp_yaw_ii*sp_pitch_d - 0.01056*sp_pitch_ii*sp_yaw_d - 0.0042*y[13]*y[8] + 0.05538*y[13]*sp_yaw_i + 0.0105*sp_yaw_ii*y[10] - 5.874*y[9] + 0.0014*sp_yaw_ii*y[7] - 0.01056*y[13]*yp[2] - 4.2e-5*y[14]*yp[1] + 5.874*sp_roll_i + 0.0042*sp_pitch_ii*y[8] - 0.0105*sp_yaw_ii*sp_pitch_i - 0.05538*sp_pitch_ii*sp_yaw_i + 0.01056*sp_pitch_ii*yp[2] + 4.2e-5*y[14]*sp_pitch_d - 0.05538*y[13]*y[11] - 0.0105*y[14]*y[10] - 0.0014*y[14]*y[7] + 0.05538*sp_pitch_ii*y[11] + 0.0105*y[14]*sp_pitch_i + 0.0001056*yp[1]*sp_yaw_d - 0.0001056*sp_yaw_d*sp_pitch_d - 4.2e-5*yp[1]*y[8] + 0.0005538*yp[1]*sp_yaw_i + 0.02641*sp_yaw_d*y[10] + 0.003521*sp_yaw_d*y[7] - 0.0001056*yp[1]*yp[2] + 4.2e-5*sp_pitch_d*y[8] - 0.02641*sp_yaw_d*sp_pitch_i - 0.0005538*sp_pitch_d*sp_yaw_i - 0.03018*y[5]*y[13] - 0.0105*y[10]*y[8] + 0.0001056*sp_pitch_d*yp[2] - 0.0014*y[7]*y[8] + 0.004221*y[1]*sp_yaw_ii + 0.1825*yaw_sp*y[13] + 0.1384*y[10]*sp_yaw_i + 0.0007*y[4]*sp_yaw_ii + 0.01846*y[7]*sp_yaw_i - 0.0005538*yp[1]*y[11] - 0.02641*y[10]*yp[2] - 0.003521*y[7]*yp[2] + 0.03018*y[5]*sp_pitch_ii + 0.0105*y[8]*sp_pitch_i - 0.0002515*thrust*y[12] - 0.1825*yaw_sp*sp_pitch_ii - 0.004221*pitch_sp*sp_yaw_ii - 0.1384*sp_yaw_i*sp_pitch_i + 0.0005538*sp_pitch_d*y[11] + 0.02641*sp_pitch_i*yp[2] + 0.0002515*thrust*sp_roll_ii - 0.004221*y[1]*y[14] - 0.1825*y[2]*y[13] - 0.1384*y[10]*y[11] - 0.0007*y[4]*y[14] - 0.01846*y[7]*y[11] + 0.004221*pitch_sp*y[14] + 0.1825*y[2]*sp_pitch_ii + 0.1384*sp_pitch_i*y[11];
              //auto k3_q = 17.44*pitch_sp - 2.892*y[4] - 17.44*y[1] + 0.001752*y[3]*sp_yaw_d - 0.0003003*y[5]*yp[0] + 0.0003003*y[5]*sp_roll_d + 0.01056*y[0]*sp_yaw_d + 0.001816*yaw_sp*yp[0] - 0.001816*yaw_sp*sp_roll_d - 0.01056*roll_sp*sp_yaw_d - 0.0006964*y[3]*y[8] - 0.01001*y[5]*y[6] + 0.009183*y[3]*sp_yaw_i + 0.06054*yaw_sp*y[6] - 0.001752*y[3]*yp[2] - 0.07507*y[5]*y[9] - 0.0042*y[0]*y[8] + 0.07507*y[5]*sp_roll_i + 0.0042*roll_sp*y[8] + 0.05537*y[0]*sp_yaw_i + 0.454*yaw_sp*y[9] - 0.01056*y[0]*yp[2] - 0.001816*y[2]*yp[0] - 0.454*yaw_sp*sp_roll_i - 0.05537*roll_sp*sp_yaw_i + 0.01056*roll_sp*yp[2] + 0.001816*y[2]*sp_roll_d - 0.009183*y[3]*y[11] - 0.06054*y[2]*y[6] - 0.05537*y[0]*y[11] - 0.454*y[2]*y[9] + 0.05537*roll_sp*y[11] + 0.454*y[2]*sp_roll_i + 0.7569*y[3]*y[5] + 0.03027*y[3]*yaw_sp - 0.03018*y[5]*y[0] + 0.03018*y[5]*roll_sp + 0.1825*y[0]*yaw_sp - 0.1825*yaw_sp*roll_sp - 0.03027*y[3]*y[2] - 0.1825*y[0]*y[2] + 0.1825*roll_sp*y[2] - 17.35*y[13] + 17.35*sp_pitch_ii + 0.004179*y[12]*sp_yaw_ii - 0.004179*sp_yaw_ii*sp_roll_ii - 0.004179*y[12]*y[14] + 0.004179*sp_roll_ii*y[14] - 0.1735*yp[1] + 0.1735*sp_pitch_d + 0.01051*y[12]*sp_yaw_d + 4.179e-5*sp_yaw_ii*yp[0] - 4.179e-5*sp_yaw_ii*sp_roll_d - 0.01051*sp_roll_ii*sp_yaw_d - 43.38*y[10] - 5.784*y[7] + 0.001393*sp_yaw_ii*y[6] - 0.004179*y[12]*y[8] + 43.38*sp_pitch_i + 0.004179*sp_roll_ii*y[8] + 0.0551*y[12]*sp_yaw_i + 0.01045*sp_yaw_ii*y[9] - 0.01051*y[12]*yp[2] - 4.179e-5*y[14]*yp[0] - 0.01045*sp_yaw_ii*sp_roll_i - 0.0551*sp_roll_ii*sp_yaw_i + 0.01051*sp_roll_ii*yp[2] + 4.179e-5*y[14]*sp_roll_d - 0.001393*y[14]*y[6] - 0.0551*y[12]*y[11] - 0.01045*y[14]*y[9] + 0.0551*sp_roll_ii*y[11] + 0.01045*y[14]*sp_roll_i + 0.0001051*yp[0]*sp_yaw_d - 0.0001051*sp_yaw_d*sp_roll_d + 0.003503*sp_yaw_d*y[6] - 4.179e-5*yp[0]*y[8] + 4.179e-5*sp_roll_d*y[8] + 0.000551*yp[0]*sp_yaw_i + 0.02627*sp_yaw_d*y[9] - 0.0001051*yp[0]*yp[2] - 0.02627*sp_yaw_d*sp_roll_i - 0.000551*sp_roll_d*sp_yaw_i - 0.001393*y[6]*y[8] + 0.0001051*sp_roll_d*yp[2] + 0.0006964*y[3]*sp_yaw_ii + 0.01837*y[6]*sp_yaw_i - 0.003503*y[6]*yp[2] - 0.03003*y[5]*y[12] - 0.01045*y[8]*y[9] + 0.03003*y[5]*sp_roll_ii + 0.01045*y[8]*sp_roll_i + 0.0042*y[0]*sp_yaw_ii + 0.1816*yaw_sp*y[12] + 0.1377*y[9]*sp_yaw_i - 0.000551*yp[0]*y[11] - 0.02627*y[9]*yp[2] - 0.1816*yaw_sp*sp_roll_ii - 0.0042*roll_sp*sp_yaw_ii - 0.1377*sp_yaw_i*sp_roll_i + 0.000551*sp_roll_d*y[11] + 0.02627*sp_roll_i*yp[2] - 0.0006964*y[3]*y[14] - 0.01837*y[6]*y[11] - 0.0042*y[0]*y[14] - 0.1816*y[2]*y[12] - 0.1377*y[9]*y[11] + 0.0042*roll_sp*y[14] + 0.1816*y[2]*sp_roll_ii + 0.1377*sp_roll_i*y[11];
              auto k3_q = 2.349*pitch_sp - 0.3896*y[4] - 2.349*y[1] + 0.001752*y[3]*sp_yaw_d - 0.0003003*y[5]*yp[0] - 2.502e-6*thrust*yp[1] + 0.0003003*y[5]*sp_roll_d + 0.01056*y[0]*sp_yaw_d + 0.001816*yaw_sp*yp[0] - 0.001816*yaw_sp*sp_roll_d - 0.01056*roll_sp*sp_yaw_d - 0.0006964*y[3]*y[8] - 0.01001*y[5]*y[6] + 2.502e-6*thrust*sp_pitch_d + 0.009183*y[3]*sp_yaw_i + 0.06054*yaw_sp*y[6] - 0.001752*y[3]*yp[2] - 0.07507*y[5]*y[9] - 0.0042*y[0]*y[8] - 0.0006255*thrust*y[10] - 8.341e-5*thrust*y[7] + 0.07507*y[5]*sp_roll_i + 0.0042*roll_sp*y[8] + 0.05537*y[0]*sp_yaw_i + 0.454*yaw_sp*y[9] - 0.01056*y[0]*yp[2] - 0.001816*y[2]*yp[0] - 0.454*yaw_sp*sp_roll_i - 0.05537*roll_sp*sp_yaw_i + 0.01056*roll_sp*yp[2] + 0.001816*y[2]*sp_roll_d + 0.0006255*thrust*sp_pitch_i - 0.009183*y[3]*y[11] - 0.06054*y[2]*y[6] - 0.05537*y[0]*y[11] - 0.454*y[2]*y[9] + 0.05537*roll_sp*y[11] + 0.454*y[2]*sp_roll_i + 0.7569*y[3]*y[5] + 0.03027*y[3]*yaw_sp - 0.03018*y[5]*y[0] - 0.0002515*y[1]*thrust - 4.17e-5*y[4]*thrust + 0.03018*y[5]*roll_sp + 0.1825*y[0]*yaw_sp - 0.1825*yaw_sp*roll_sp + 0.0002515*pitch_sp*thrust - 0.03027*y[3]*y[2] - 0.1825*y[0]*y[2] + 0.1825*roll_sp*y[2] - 2.338*y[13] + 2.338*sp_pitch_ii + 0.004179*y[12]*sp_yaw_ii - 0.004179*sp_yaw_ii*sp_roll_ii - 0.004179*y[12]*y[14] + 0.004179*sp_roll_ii*y[14] - 0.02338*yp[1] + 0.02338*sp_pitch_d + 0.01051*y[12]*sp_yaw_d + 4.179e-5*sp_yaw_ii*yp[0] - 4.179e-5*sp_yaw_ii*sp_roll_d - 0.01051*sp_roll_ii*sp_yaw_d - 5.844*y[10] - 0.7792*y[7] + 0.001393*sp_yaw_ii*y[6] - 0.004179*y[12]*y[8] + 5.844*sp_pitch_i + 0.004179*sp_roll_ii*y[8] + 0.0551*y[12]*sp_yaw_i + 0.01045*sp_yaw_ii*y[9] - 0.01051*y[12]*yp[2] - 4.179e-5*y[14]*yp[0] - 0.01045*sp_yaw_ii*sp_roll_i - 0.0551*sp_roll_ii*sp_yaw_i + 0.01051*sp_roll_ii*yp[2] + 4.179e-5*y[14]*sp_roll_d - 0.001393*y[14]*y[6] - 0.0551*y[12]*y[11] - 0.01045*y[14]*y[9] + 0.0551*sp_roll_ii*y[11] + 0.01045*y[14]*sp_roll_i + 0.0001051*yp[0]*sp_yaw_d - 0.0001051*sp_yaw_d*sp_roll_d + 0.003503*sp_yaw_d*y[6] - 4.179e-5*yp[0]*y[8] + 4.179e-5*sp_roll_d*y[8] + 0.000551*yp[0]*sp_yaw_i + 0.02627*sp_yaw_d*y[9] - 0.0001051*yp[0]*yp[2] - 0.02627*sp_yaw_d*sp_roll_i - 0.000551*sp_roll_d*sp_yaw_i - 0.001393*y[6]*y[8] + 0.0001051*sp_roll_d*yp[2] + 0.0006964*y[3]*sp_yaw_ii + 0.01837*y[6]*sp_yaw_i - 0.003503*y[6]*yp[2] - 0.03003*y[5]*y[12] - 0.01045*y[8]*y[9] - 0.0002502*thrust*y[13] + 0.03003*y[5]*sp_roll_ii + 0.01045*y[8]*sp_roll_i + 0.0042*y[0]*sp_yaw_ii + 0.1816*yaw_sp*y[12] + 0.1377*y[9]*sp_yaw_i - 0.000551*yp[0]*y[11] - 0.02627*y[9]*yp[2] - 0.1816*yaw_sp*sp_roll_ii - 0.0042*roll_sp*sp_yaw_ii - 0.1377*sp_yaw_i*sp_roll_i + 0.000551*sp_roll_d*y[11] + 0.02627*sp_roll_i*yp[2] + 0.0002502*thrust*sp_pitch_ii - 0.0006964*y[3]*y[14] - 0.01837*y[6]*y[11] - 0.0042*y[0]*y[14] - 0.1816*y[2]*y[12] - 0.1377*y[9]*y[11] + 0.0042*roll_sp*y[14] + 0.1816*y[2]*sp_roll_ii + 0.1377*sp_roll_i*y[11];
              
              auto k1_r = -5.427e-8; // p_dot * q_dot coefficient
              auto k2_r = 3.273e-5*pitch_sp - 5.427e-6*y[4] - 3.273e-5*y[1] - 3.256e-5*y[13] + 3.256e-5*sp_pitch_ii - 3.256e-7*yp[1] + 3.256e-7*sp_pitch_d - 8.141e-5*y[10] - 1.085e-5*y[7] + 8.141e-5*sp_pitch_i; // p_dot coefficient
              auto k3_r = 3.273e-5*roll_sp - 3.273e-5*y[0] - 5.427e-6*y[3] - 3.256e-5*y[12] + 3.256e-5*sp_roll_ii - 3.256e-7*yp[0] + 3.256e-7*sp_roll_d - 1.085e-5*y[6] - 8.141e-5*y[9] + 8.141e-5*sp_roll_i; // q_dot coefficient
              //auto k4_r = 1.748*yaw_sp - 0.289*y[5] - 1.748*y[2] - 3.256e-5*y[3]*yp[1] - 0.0001964*y[1]*yp[0] - 0.0001964*y[0]*yp[1] - 3.256e-5*y[4]*yp[0] + 3.256e-5*y[3]*sp_pitch_d + 0.0001964*y[1]*sp_roll_d + 0.0001964*roll_sp*yp[1] + 3.256e-5*y[4]*sp_roll_d - 0.008141*y[3]*y[10] - 0.006545*y[1]*y[6] + 0.0001964*y[0]*sp_pitch_d + 0.0001964*pitch_sp*yp[0] - 0.001085*y[3]*y[7] - 0.001085*y[4]*y[6] - 0.0001964*roll_sp*sp_pitch_d - 0.0001964*pitch_sp*sp_roll_d - 0.04909*y[1]*y[9] - 0.04909*y[0]*y[10] - 0.008141*y[4]*y[9] - 0.006545*y[0]*y[7] + 0.008141*y[3]*sp_pitch_i + 0.006545*pitch_sp*y[6] + 0.04909*y[1]*sp_roll_i + 0.04909*roll_sp*y[10] + 0.008141*y[4]*sp_roll_i + 0.006545*roll_sp*y[7] + 0.04909*y[0]*sp_pitch_i + 0.04909*pitch_sp*y[9] - 0.04909*roll_sp*sp_pitch_i - 0.04909*pitch_sp*sp_roll_i - 0.003273*y[3]*y[1] - 0.00341*y[3]*y[4] - 0.01973*y[1]*y[0] - 0.003273*y[4]*y[0] + 0.003273*y[3]*pitch_sp + 0.01973*y[1]*roll_sp + 0.003273*y[4]*roll_sp + 0.01973*y[0]*pitch_sp - 0.01973*roll_sp*pitch_sp - 0.01954*y[13]*y[12] + 0.01954*y[13]*sp_roll_ii + 0.04022*sp_yaw_ii + 0.01954*y[12]*sp_pitch_ii - 0.01954*sp_roll_ii*sp_pitch_ii - 0.04022*y[14] - 0.0001954*y[13]*yp[0] - 0.0001954*y[12]*yp[1] + 0.0001954*y[13]*sp_roll_d + 0.0001954*sp_roll_ii*yp[1] + 0.1012*sp_yaw_d - 0.006513*y[13]*y[6] + 0.0001954*y[12]*sp_pitch_d + 0.0001954*sp_pitch_ii*yp[0] - 0.0001954*sp_roll_ii*sp_pitch_d - 0.0001954*sp_pitch_ii*sp_roll_d - 0.04884*y[13]*y[9] - 0.04884*y[12]*y[10] - 0.006513*y[12]*y[7] + 0.006513*sp_pitch_ii*y[6] - 0.04022*y[8] + 0.04884*y[13]*sp_roll_i + 0.04884*sp_roll_ii*y[10] + 0.006513*sp_roll_ii*y[7] + 0.5304*sp_yaw_i + 0.04884*y[12]*sp_pitch_i + 0.04884*sp_pitch_ii*y[9] - 0.1012*yp[2] - 0.04884*sp_roll_ii*sp_pitch_i - 0.04884*sp_pitch_ii*sp_roll_i - 0.5304*y[11] - 1.954e-6*yp[1]*yp[0] + 1.954e-6*yp[1]*sp_roll_d - 6.513e-5*yp[1]*y[6] + 1.954e-6*yp[0]*sp_pitch_d - 1.954e-6*sp_roll_d*sp_pitch_d - 0.0004884*yp[1]*y[9] - 0.0004884*yp[0]*y[10] - 6.513e-5*yp[0]*y[7] + 6.513e-5*sp_pitch_d*y[6] + 0.0004884*yp[1]*sp_roll_i + 0.0004884*sp_roll_d*y[10] + 6.513e-5*sp_roll_d*y[7] - 0.003256*y[3]*y[13] - 0.01628*y[6]*y[10] + 0.0004884*yp[0]*sp_pitch_i + 0.0004884*sp_pitch_d*y[9] - 0.002171*y[6]*y[7] - 0.0004884*sp_roll_d*sp_pitch_i - 0.0004884*sp_pitch_d*sp_roll_i - 0.01964*y[1]*y[12] - 0.01964*y[0]*y[13] - 0.1221*y[10]*y[9] - 0.003256*y[4]*y[12] - 0.01628*y[7]*y[9] + 0.003256*y[3]*sp_pitch_ii + 0.01628*y[6]*sp_pitch_i + 0.01964*y[1]*sp_roll_ii + 0.01964*roll_sp*y[13] + 0.1221*y[10]*sp_roll_i + 0.003256*y[4]*sp_roll_ii + 0.01628*y[7]*sp_roll_i + 0.01964*y[0]*sp_pitch_ii + 0.01964*pitch_sp*y[12] + 0.1221*y[9]*sp_pitch_i - 0.01964*roll_sp*sp_pitch_ii - 0.01964*pitch_sp*sp_roll_ii - 0.1221*sp_roll_i*sp_pitch_i;
              auto k4_r = 0.2355*yaw_sp - 0.03894*y[5] - 0.2355*y[2] - 3.256e-5*y[3]*yp[1] - 0.0001964*y[1]*yp[0] - 0.0001964*y[0]*yp[1] - 3.256e-5*y[4]*yp[0] + 3.256e-5*y[3]*sp_pitch_d + 0.0001964*y[1]*sp_roll_d + 0.0001964*roll_sp*yp[1] + 3.256e-5*y[4]*sp_roll_d - 0.008141*y[3]*y[10] - 0.006545*y[1]*y[6] + 0.0001964*y[0]*sp_pitch_d + 0.0001964*pitch_sp*yp[0] - 0.001085*y[3]*y[7] - 0.001085*y[4]*y[6] - 0.0001964*roll_sp*sp_pitch_d - 0.0001964*pitch_sp*sp_roll_d + 1.459e-6*thrust*sp_yaw_d - 0.04909*y[1]*y[9] - 0.04909*y[0]*y[10] - 0.008141*y[4]*y[9] - 0.006545*y[0]*y[7] + 0.008141*y[3]*sp_pitch_i + 0.006545*pitch_sp*y[6] + 0.04909*y[1]*sp_roll_i + 0.04909*roll_sp*y[10] + 0.008141*y[4]*sp_roll_i + 0.006545*roll_sp*y[7] + 0.04909*y[0]*sp_pitch_i + 0.04909*pitch_sp*y[9] - 5.8e-7*thrust*y[8] - 0.04909*roll_sp*sp_pitch_i - 0.04909*pitch_sp*sp_roll_i + 7.648e-6*thrust*sp_yaw_i - 1.459e-6*thrust*yp[2] - 7.648e-6*thrust*y[11] - 0.003273*y[3]*y[1] - 0.00341*y[3]*y[4] - 0.01973*y[1]*y[0] - 0.003273*y[4]*y[0] + 0.003273*y[3]*pitch_sp + 0.01973*y[1]*roll_sp + 0.003273*y[4]*roll_sp + 0.01973*y[0]*pitch_sp - 4.168e-6*y[5]*thrust - 0.01973*roll_sp*pitch_sp + 2.521e-5*yaw_sp*thrust - 2.521e-5*thrust*y[2] - 0.01954*y[13]*y[12] + 0.01954*y[13]*sp_roll_ii + 0.005419*sp_yaw_ii + 0.01954*y[12]*sp_pitch_ii - 0.01954*sp_roll_ii*sp_pitch_ii - 0.005419*y[14] - 0.0001954*y[13]*yp[0] - 0.0001954*y[12]*yp[1] + 0.0001954*y[13]*sp_roll_d + 0.0001954*sp_roll_ii*yp[1] + 0.01363*sp_yaw_d - 0.006513*y[13]*y[6] + 0.0001954*y[12]*sp_pitch_d + 0.0001954*sp_pitch_ii*yp[0] - 0.0001954*sp_roll_ii*sp_pitch_d - 0.0001954*sp_pitch_ii*sp_roll_d - 0.04884*y[13]*y[9] - 0.04884*y[12]*y[10] - 0.006513*y[12]*y[7] + 0.006513*sp_pitch_ii*y[6] - 0.005419*y[8] + 0.04884*y[13]*sp_roll_i + 0.04884*sp_roll_ii*y[10] + 0.006513*sp_roll_ii*y[7] + 0.07145*sp_yaw_i + 0.04884*y[12]*sp_pitch_i + 0.04884*sp_pitch_ii*y[9] - 0.01363*yp[2] - 0.04884*sp_roll_ii*sp_pitch_i - 0.04884*sp_pitch_ii*sp_roll_i - 0.07145*y[11] - 1.954e-6*yp[1]*yp[0] + 1.954e-6*yp[1]*sp_roll_d - 6.513e-5*yp[1]*y[6] + 1.954e-6*yp[0]*sp_pitch_d - 1.954e-6*sp_roll_d*sp_pitch_d - 0.0004884*yp[1]*y[9] - 0.0004884*yp[0]*y[10] - 6.513e-5*yp[0]*y[7] + 6.513e-5*sp_pitch_d*y[6] + 0.0004884*yp[1]*sp_roll_i + 0.0004884*sp_roll_d*y[10] + 6.513e-5*sp_roll_d*y[7] - 0.003256*y[3]*y[13] - 0.01628*y[6]*y[10] + 0.0004884*yp[0]*sp_pitch_i + 0.0004884*sp_pitch_d*y[9] - 0.002171*y[6]*y[7] - 0.0004884*sp_roll_d*sp_pitch_i - 0.0004884*sp_pitch_d*sp_roll_i - 0.01964*y[1]*y[12] - 0.01964*y[0]*y[13] - 0.1221*y[10]*y[9] - 0.003256*y[4]*y[12] - 0.01628*y[7]*y[9] + 0.003256*y[3]*sp_pitch_ii + 0.01628*y[6]*sp_pitch_i + 0.01964*y[1]*sp_roll_ii + 0.01964*roll_sp*y[13] + 0.1221*y[10]*sp_roll_i + 0.003256*y[4]*sp_roll_ii + 0.01628*y[7]*sp_roll_i + 0.01964*y[0]*sp_pitch_ii + 0.01964*pitch_sp*y[12] + 0.1221*y[9]*sp_pitch_i - 0.01964*roll_sp*sp_pitch_ii - 0.01964*pitch_sp*sp_roll_ii - 0.1221*sp_roll_i*sp_pitch_i + 5.8e-7*thrust*sp_yaw_ii - 5.8e-7*thrust*y[14];
              
              yp[3] = (k3_p * (1.0 - k2_q) + k2_p*k3_q)/denom;
              yp[4] = (k3_q * (1.0 - k1_p) + k1_q*k3_p)/denom;
              yp[5] = k1_r * yp[3] * yp[4] + k2_r * yp[3] + k3_r * yp[4] + k4_r;
          } else if(syschoice == 11)
          {
              /* Position control using symbolic calculus in matlab */
              /* The only uncertaincy are from the initial states */
              /* All the physical parameters are certain */
              static const double x_sp = 0.0;
              static const double y_sp = 0.0;
              static const double z_sp = 1.0;
              static const double yaw_sp = 0.0;
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = 1.0;//cos(y[2]*M_PI/180.0);
              auto sinYaw = 0.0;//sin(y[2]*M_PI/180.0);
              
              auto tanPitch = tan(y[1]);
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = 0;//(y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // Translation kinematic x, y , z
              yp[15] = y[20]*(sinRoll*sinYaw + cosRoll*cosYaw*sinPitch) - y[19]*(cosRoll*sinYaw - cosYaw*sinPitch*sinRoll) + cosPitch*cosYaw*y[18];
              yp[16] = y[19]*(cosRoll*cosYaw + sinPitch*sinRoll*sinYaw) - y[20]*(cosYaw*sinRoll - cosRoll*sinPitch*sinYaw) + cosPitch*sinYaw*y[18];
              yp[17] = cosPitch*cosRoll*y[20] - sinPitch*y[18] + cosPitch*sinRoll*y[19];
              
              // Translation kinematic u , v , w
              yp[18] = y[5]*y[19] - y[4]*y[20] + g*sinPitch;
              yp[19] = y[3]*y[20] - y[5]*y[18] - g*cosPitch*sinRoll;
              
              // integral of error in x , y and z
              auto e_X = x_sp - y[15];
              auto e_Y = y_sp - y[16];
              auto e_Z = z_sp - y[17];
              yp[21] = 2*e_X - y[18];
              yp[22] = 2*e_Y - y[19];
              yp[23] = 2*e_Z - y[20];
              
              // integral of speed u , v and w
              //yp[24] = y[18];
              //yp[25] = y[19];
              //yp[26] = y[20];
              
              auto temp1 = 50*e_X - 25*y[18] + y[21];
              auto temp2 = 50*e_Y - 25*y[19] + y[22];
              auto roll_sp = (sinYaw*temp1 - cosYaw*temp2);
              auto pitch_sp = -(- cosYaw*temp1 - sinYaw*temp2);
              // auto roll_sp = temp2;//sinYaw*temp1 - cosYaw*temp2;
              // auto pitch_sp = -temp1;
              
              // integral of error in roll , pitch and yaw
              yp[6] = roll_sp - y[0];
              yp[7] = pitch_sp - y[1];
              yp[8] = yaw_sp - y[2];
              
              // Double integrale of error in roll , pitch and yaw
              yp[9] = y[6];
              yp[10] = y[7];
              yp[11] = y[8];
              
              // integral of p , r and q
              yp[12] = y[3];
              yp[13] = y[5];
              yp[14] = y[4];
              
              auto diff_eX = - yp[15];
              auto diff_eY = - yp[16];
              auto diff_eZ = - yp[17];
              
              
              auto y1 = -yp[2];
              auto p1 = -((sinYaw*temp1 - cosYaw*temp2)*yp[2] - sinYaw*(2*e_Y - y[19] + 50*diff_eY - 25*yp[19]) - cosYaw*(2*e_X - y[18] + 50*diff_eX - 25*yp[18])) - yp[1];
              auto r1 = (sinYaw*(2*e_X - y[18] + 50*diff_eX - 25*yp[18]) - cosYaw*(2*e_Y - y[19] + 50*diff_eY - 25*yp[19]) + (cosYaw*temp1 + sinYaw*temp2)*yp[2]) - yp[0]; // diff yaw_sp - yaw
              // auto p1 = -(yp[21] + 50*diff_eX - 25*yp[18]) - yp[1];
              // auto r1 = (yp[22] + 50*diff_eY - 25*yp[19])  - yp[0]; // diff yaw_sp - yaw
              
              auto thrust = 50000*e_Z - 25000*y[20] + 15000*y[23] + 36000;
              
              auto k1_p = - 4.202e-7*thrust - 0.003926;
              auto k2_p = 0.000305*yp[8] - 5.043e-5*y[5] + 1.765e-5*y1 + 7.018e-6*y[11] + 9.253e-5*y[8] - 7.018e-6*y[13];
              
              auto k1_q = 0.0003055*yp[8] - 5.051e-5*y[5] + 1.768e-5*y1 + 7.029e-6*y[11] + 9.269e-5*y[8] - 7.029e-6*y[13];
              auto k2_q = - 4.209e-7*thrust - 0.003932;
              
              auto denom  = (1.0 - k1_p) * (1.0 - k2_q) - k1_q*k2_p;
              auto k3_p = 2.367*yp[6] - 0.3926*y[3] + 0.02355*r1 - 0.02647*y1*y[7] + 4.211e-5*p1*y[13] - 0.0007018*y[4]*y[13] - 0.01009*y[5]*y[14] - 8.404e-5*thrust*y[12] + 0.00353*y1*y[14] - 0.1839*yp[7]*yp[8] - 0.00183*yp[8]*p1 + 0.0305*yp[8]*y[4] + 0.03041*yp[7]*y[5] + 0.0002534*yp[6]*thrust - 0.01064*yp[7]*y1 + 0.0003026*p1*y[5] - 1.045*y[4]*y[5] - 4.202e-5*y[3]*thrust + 2.521e-6*r1*thrust - 0.0001059*p1*y1 + 0.001765*y[4]*y1 - 0.004211*y[10]*y[11] + 2.355*y[9] - 0.05552*y[10]*y[8] - 0.01053*y[11]*y[7] + 0.001404*y[11]*y[14] + 5.889*y[6] + 0.004211*y[10]*y[13] - 0.7852*y[12] - 0.004232*yp[7]*y[11] - 0.183*yp[8]*y[10] - 0.1388*y[7]*y[8] - 4.211e-5*p1*y[11] + 0.0007018*y[4]*y[11] + 0.01851*y[8]*y[14] + 0.03026*y[5]*y[10] + 0.01053*y[7]*y[13] + 0.0002521*thrust*y[9] - 0.01059*y1*y[10] - 0.001404*y[14]*y[13] - 0.0558*yp[7]*y[8] - 0.4575*yp[8]*y[7] - 0.0005552*p1*y[8] + 0.061*yp[8]*y[14] + 0.009253*y[4]*y[8] + 0.004232*yp[7]*y[13] + 0.07564*y[5]*y[7] + 0.0006303*thrust*y[6];
              auto k3_q = 2.371*yp[7] + 0.02359*p1 - 0.3932*y[4] - 0.02652*y1*y[6] - 0.0007029*y[3]*y[13] - 0.0101*y[5]*y[12] + 4.218e-5*r1*y[13] - 8.418e-5*thrust*y[14] + 0.003536*y1*y[12] - 0.1842*yp[6]*yp[8] + 0.03055*yp[8]*y[3] + 0.03046*yp[6]*y[5] - 0.001833*yp[8]*r1 + 0.0002538*yp[7]*thrust - 0.01066*yp[6]*y1 + 1.035*y[3]*y[5] + 0.0003031*y[5]*r1 + 2.526e-6*p1*thrust - 4.209e-5*y[4]*thrust + 0.001768*y[3]*y1 - 0.0001061*r1*y1 - 0.004218*y[9]*y[11] + 2.359*y[10] - 0.05561*y[9]*y[8] - 0.01054*y[11]*y[6] + 0.001406*y[11]*y[12] + 5.898*y[7] + 0.004218*y[9]*y[13] - 0.7865*y[14] - 0.004239*yp[6]*y[11] - 0.1833*yp[8]*y[9] - 0.139*y[6]*y[8] + 0.0007029*y[3]*y[11] + 0.01854*y[8]*y[12] + 0.03031*y[5]*y[9] + 0.01054*y[6]*y[13] - 4.218e-5*r1*y[11] + 0.0002526*thrust*y[10] - 0.01061*y1*y[9] - 0.001406*y[12]*y[13] - 0.05589*yp[6]*y[8] - 0.4583*yp[8]*y[6] + 0.0611*yp[8]*y[12] + 0.009269*y[3]*y[8] + 0.004239*yp[6]*y[13] + 0.07577*y[5]*y[6] - 0.0005561*r1*y[8] + 0.0006314*thrust*y[7];
              
              //auto k1_r = -5.466e-8; // p_dot * q_dot coefficient
              //auto k2_r = 3.296e-5*yp[7] + 3.28e-7*p1 - 5.466e-6*y[4] + 3.28e-5*y[10] + 8.199e-5*y[7] - 1.093e-5*y[14]; // p_dot coefficient
              //auto k3_r = 3.296e-5*yp[6] - 5.466e-6*y[3] + 3.28e-7*r1 + 3.28e-5*y[9] + 8.199e-5*y[6] - 1.093e-5*y[12]; // q_dot coefficient
              auto k4_r = 0.2372*yp[8] - 0.03922*y[5] + 0.01373*y1 + 6.559e-5*p1*y[12] - 0.001093*y[3]*y[14] - 0.001093*y[4]*y[12] + 6.559e-5*r1*y[14] - 5.842e-7*thrust*y[13] - 0.01988*yp[7]*yp[6] + 0.003296*yp[7]*y[3] - 0.0001978*yp[6]*p1 + 0.003296*yp[6]*y[4] - 0.0001978*yp[7]*r1 + 2.539e-5*yp[8]*thrust + 3.28e-5*y[3]*p1 - 0.0005466*y[3]*y[4] - 0.015*y[3]*y[5] - 1.968e-6*p1*r1 + 3.28e-5*y[4]*r1 - 4.198e-6*y[5]*thrust + 1.469e-6*thrust*y1 - 0.01968*y[10]*y[9] + 0.005458*y[11] - 0.04919*y[10]*y[6] - 0.04919*y[9]*y[7] + 0.006559*y[10]*y[12] + 0.006559*y[9]*y[14] + 0.07197*y[8] - 0.005458*y[13] - 0.01978*yp[7]*y[9] - 0.01978*yp[6]*y[10] - 0.123*y[7]*y[6] + 0.00328*y[3]*y[10] + 0.0164*y[7]*y[12] - 0.0001968*p1*y[9] + 0.00328*y[4]*y[9] + 0.0164*y[6]*y[14] - 0.0001968*r1*y[10] + 5.842e-7*thrust*y[11] - 0.002186*y[12]*y[14] - 0.04944*yp[7]*y[6] - 0.04944*yp[6]*y[7] + 0.006592*yp[7]*y[12] + 0.008199*y[3]*y[7] - 0.0004919*p1*y[6] + 0.006592*yp[6]*y[14] + 0.008199*y[4]*y[6] - 0.0004919*r1*y[7] + 7.703e-6*thrust*y[8]; // other
              
              // cout << "DEnom : " << getAAF(thrust).convert_int() << endl;
              // cout << "DEnom : " << thrust << endl;
              yp[3] = (k3_p * (1.0 - k2_q) + k2_p*k3_q)/denom;
              yp[4] = (k3_q * (1.0 - k1_p) + k1_q*k3_p)/denom;
              yp[5] = /*k1_r * yp[3] * yp[4] + k2_r * yp[3] + k3_r * yp[4]*/ k4_r;
              
              //Z coordinate
              auto F = 1.596e-6*thrust + 2.669e-8*y[4]*yp[4] + 5.339e-6*y[3]*y[12] + 5.339e-6*y[4]*y[14] + 3.423e-7*y[5]*y[13] + 4.805e-9*p1*p1 + 4.805e-9*r1*r1 + 1.507e-7*y1*y1 + 0.0003003*y[7]*y[7] + 0.0003003*y[6]*y[6] + 4.142e-6*y[8]*y[8] + 1.335e-10*yp[3]*yp[3] - 1.61e-5*yp[6]*y[3] + 1.335e-10*yp[4]*yp[4] - 1.61e-5*yp[7]*y[4] - 1.488e-5*yp[8]*y[5] + 5.339e-6*y[12]*y[12] + 5.339e-6*y[14]*y[14] + 2.382e-8*y[13]*y[13] + 4.853e-5*yp[7]*yp[7] + 4.853e-5*yp[6]*yp[6] + 4.5e-5*yp[8]*yp[8] + 1.335e-6*y[3]*y[3] + 1.335e-6*y[4]*y[4] + 1.23e-6*y[5]*y[5] + 8.542e-11*thrust*thrust + 9.609e-7*y[10]*p1 + 9.609e-7*y[9]*r1 + 1.198e-7*y[11]*y1 + 0.0002402*y[10]*y[7] + 0.0002402*y[9]*y[6] + 6.282e-7*y[11]*y[8] - 1.602e-7*y[9]*yp[3] - 1.602e-7*y[10]*yp[4] - 3.203e-5*y[9]*y[12] - 3.203e-5*y[10]*y[14] - 4.764e-8*y[11]*y[13] + 2.402e-6*p1*y[7] + 2.402e-6*r1*y[6] + 1.58e-6*y1*y[8] + 9.657e-5*yp[7]*y[10] + 9.657e-5*yp[6]*y[9] + 2.071e-6*yp[8]*y[11] - 1.602e-9*r1*yp[3] - 1.602e-9*p1*yp[4] - 3.203e-7*r1*y[12] - 4.004e-7*y[6]*yp[3] - 3.203e-7*p1*y[14] - 4.004e-7*y[7]*yp[4] - 1.198e-7*y1*y[13] - 1.602e-5*y[3]*y[9] - 8.008e-5*y[6]*y[12] - 1.602e-5*y[4]*y[10] - 8.008e-5*y[7]*y[14] - 3.423e-7*y[5]*y[11] - 6.282e-7*y[8]*y[13] + 5.339e-8*yp[3]*y[12] + 5.339e-8*yp[4]*y[14] + 9.657e-7*yp[7]*p1 + 9.657e-7*yp[6]*r1 + 5.208e-6*yp[8]*y1 + 0.0002414*yp[7]*y[7] + 0.0002414*yp[6]*y[6] + 2.73e-5*yp[8]*y[8] + 4.805e-5*y[10]*y[10] + 4.805e-5*y[9]*y[9] + 2.382e-8*y[11]*y[11] - 1.61e-7*yp[6]*yp[3] - 1.602e-7*y[3]*r1 - 1.61e-7*yp[7]*yp[4] - 1.602e-7*y[4]*p1 - 8.61e-7*y[5]*y1 - 3.219e-5*yp[6]*y[12] - 4.004e-5*y[3]*y[6] - 3.219e-5*yp[7]*y[14] - 4.004e-5*y[4]*y[7] - 2.071e-6*yp[8]*y[13] - 4.514e-6*y[5]*y[8] + 2.669e-8*y[3]*yp[3] + 0.007455;
              //auto F = 1.596e-6*thrust  + 5.339e-6*y[3]*y[12] + 5.339e-6*y[4]*y[14] + 3.423e-7*y[5]*y[13]  + 0.0003003*y[7]*y[7] + 0.0003003*y[6]*y[6] + 4.142e-6*y[8]*y[8] - 1.61e-5*yp[6]*y[3] - 1.61e-5*yp[7]*y[4] - 1.488e-5*yp[8]*y[5] + 5.339e-6*y[12]*y[12] + 5.339e-6*y[14]*y[14]  + 4.853e-5*yp[7]*yp[7] + 4.853e-5*yp[6]*yp[6] + 4.5e-5*yp[8]*yp[8] + 1.335e-6*y[3]*y[3] + 1.335e-6*y[4]*y[4] + 1.23e-6*y[5]*y[5] + 8.542e-11*thrust*thrust + 9.609e-7*y[10]*p1 + 0.0002402*y[10]*y[7] + 0.0002402*y[9]*y[6] - 3.203e-5*y[9]*y[12] - 3.203e-5*y[10]*y[14] - 4.764e-8*y[11]*y[13] + 2.402e-6*p1*y[7] + 2.402e-6*r1*y[6] + 1.58e-6*y1*y[8] + 9.657e-5*yp[7]*y[10] + 9.657e-5*yp[6]*y[9] + 2.071e-6*yp[8]*y[11] - 1.602e-9*r1*yp[3]  - 1.602e-5*y[3]*y[9] - 8.008e-5*y[6]*y[12] - 1.602e-5*y[4]*y[10] - 8.008e-5*y[7]*y[14] - 3.423e-7*y[5]*y[11] - 6.282e-7*y[8]*y[13] + 5.339e-8*yp[3]*y[12] + 5.339e-8*yp[4]*y[14] + 9.657e-7*yp[7]*p1 + 9.657e-7*yp[6]*r1 + 5.208e-6*yp[8]*y1 + 0.0002414*yp[7]*y[7] + 0.0002414*yp[6]*y[6] + 2.73e-5*yp[8]*y[8] + 4.805e-5*y[10]*y[10] + 4.805e-5*y[9]*y[9] - 3.219e-5*yp[6]*y[12] - 4.004e-5*y[3]*y[6] - 3.219e-5*yp[7]*y[14] - 4.004e-5*y[4]*y[7] - 2.071e-6*yp[8]*y[13] - 4.514e-6*y[5]*y[8] + 0.007455;
              yp[20] = y[4]*y[18] - y[3]*y[19] + F/m - g*cosPitch*cosRoll;
              
          } else if (syschoice == 12) // We test here height and rpy rate control only
          {
              /*Only uncertaincy on the intial states-> params are not*/
              /*Derivative terms are (Kd coefficient are zro (?) */
              static const double z_sp = 1.0;
              static const double q_sp = 0.0;
              static const double p_sp = 1.0;
              static const double r_sp = 0.0;
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              // Because our target are constant the integrale are val * t
              auto int_rsp = r_sp * y[14];
              auto int_psp = p_sp * y[14];
              auto int_qsp = q_sp * y[14];
              
              // You can derive derivative of q_sp , r_sp or p_sp
              // If they are not constant
              auto d_q_sp = 0;
              auto d_p_sp = 0;
              auto d_r_sp = 0;
              
              auto thrust_sp = 50000*z_sp - 50000*y[12] - 25000*y[11] - 15000*y[13] + 36000;
              // cout << getAAF(thrust_sp).convert_int() << endl;
              // cout << getAAF(z_sp - y[12]).convert_int() << endl;
              // cout << "------------" << endl;
              
              auto k1_p = - 4.191e-7*thrust_sp - 0.003916;
              auto k2_p = 5.03e-5*r_sp - 5.03e-5*y[5] - 7.0e-6*y[8] + 7.0e-6*int_rsp;
              
              auto k1_q = 5.004e-5*r_sp - 5.004e-5*y[5] - 6.964e-6*y[8] + 6.964e-6*int_rsp;
              auto k2_q = - 4.17e-7*thrust_sp - 0.003896;
              
              auto denom  = (1.0 - k1_p) * (1.0 - k2_q) - k1_q*k2_p;
              auto k3_p = 0.3916*p_sp - 0.3916*y[3] + 5.03e-5*y[5]*d_q_sp - 5.03e-5*r_sp*d_q_sp + 4.191e-7*thrust_sp*d_p_sp - 0.0007*y[4]*y[8] - 0.01006*y[5]*y[7] + 0.0007*y[4]*int_rsp + 0.0007*q_sp*y[8] + 0.01006*y[5]*int_qsp + 0.01006*r_sp*y[7] - 0.0007*q_sp*int_rsp - 0.01006*r_sp*int_qsp - 8.383e-5*thrust_sp*y[6] + 8.383e-5*thrust_sp*int_psp - 0.7657*y[4]*y[5] + 0.00503*y[4]*r_sp + 0.00503*q_sp*y[5] - 0.00503*q_sp*r_sp - 4.191e-5*y[3]*thrust_sp + 4.191e-5*p_sp*thrust_sp + 0.003916*d_p_sp - 0.7831*y[6] + 0.7831*int_psp + 7.0e-6*d_q_sp*y[8] - 7.0e-6*d_q_sp*int_rsp - 0.0014*y[7]*y[8] + 0.0014*y[7]*int_rsp + 0.0014*int_qsp*y[8] - 0.0014*int_qsp*int_rsp;
              auto k3_q = 0.3896*q_sp - 0.3896*y[4] + 5.004e-5*y[5]*d_p_sp - 5.004e-5*r_sp*d_p_sp + 4.17e-7*thrust_sp*d_q_sp - 0.0006964*y[3]*y[8] - 0.01001*y[5]*y[6] + 0.0006964*y[3]*int_rsp + 0.0006964*p_sp*y[8] + 0.01001*y[5]*int_psp + 0.01001*r_sp*y[6] - 0.0006964*p_sp*int_rsp - 0.01001*r_sp*int_psp - 8.341e-5*thrust_sp*y[7] + 8.341e-5*thrust_sp*int_qsp + 0.7569*y[3]*y[5] + 0.005004*y[3]*r_sp + 0.005004*p_sp*y[5] - 0.005004*p_sp*r_sp - 4.17e-5*y[4]*thrust_sp + 4.17e-5*q_sp*thrust_sp + 0.003896*d_q_sp - 0.7792*y[7] + 0.7792*int_qsp + 6.964e-6*d_p_sp*y[8] - 6.964e-6*d_p_sp*int_rsp - 0.001393*y[6]*y[8] + 0.001393*y[6]*int_rsp + 0.001393*int_psp*y[8] - 0.001393*int_psp*int_rsp;
              
              // auto k1_r = -5.466e-8; // p_dot * q_dot coefficient
              auto k2_r = 5.427e-6*q_sp - 5.427e-6*y[4] + 5.427e-8*d_q_sp - 1.085e-5*y[7] + 1.085e-5*int_qsp; // p_dot coefficient
              auto k3_r = 5.427e-6*p_sp - 5.427e-6*y[3] + 5.427e-8*d_p_sp - 1.085e-5*y[6] + 1.085e-5*int_psp; // q_dot coefficient
              auto k4_r = 0.03894*r_sp - 0.03894*y[5] + 5.427e-6*y[3]*d_q_sp + 5.427e-6*y[4]*d_p_sp - 5.427e-6*p_sp*d_q_sp - 5.427e-6*q_sp*d_p_sp - 0.001085*y[3]*y[7] - 0.001085*y[4]*y[6] + 0.001085*y[3]*int_qsp + 0.001085*p_sp*y[7] + 0.001085*y[4]*int_psp + 0.001085*q_sp*y[6] - 0.001085*p_sp*int_qsp - 0.001085*q_sp*int_psp - 5.8e-7*thrust_sp*y[8] + 5.8e-7*thrust_sp*int_rsp - 0.00341*y[3]*y[4] + 0.0005427*y[3]*q_sp + 0.0005427*p_sp*y[4] - 0.0005427*p_sp*q_sp - 4.168e-6*y[5]*thrust_sp + 4.168e-6*r_sp*thrust_sp - 0.005419*y[8] + 0.005419*int_rsp - 5.427e-8*d_p_sp*d_q_sp + 1.085e-5*d_p_sp*y[7] + 1.085e-5*d_q_sp*y[6] - 1.085e-5*d_p_sp*int_qsp - 1.085e-5*d_q_sp*int_psp - 0.002171*y[6]*y[7] + 0.002171*y[6]*int_qsp + 0.002171*int_psp*y[7] - 0.002171*int_psp*int_qsp;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q , r
              yp[3] = (k3_p * (1.0 - k2_q) + k2_p*k3_q)/denom;
              yp[4] = (k3_q * (1.0 - k1_p) + k1_q*k3_p)/denom;
              yp[5] = /*k1_r * yp[3] * yp[4] +*/ k2_r * yp[3] + k3_r * yp[4] + k4_r;
              
              // Integral of p , q and r
              yp[6] = y[3];
              yp[7] = y[4];
              yp[8] = y[5];
              
              // Translation kinematic u , v , w
              auto F = 1.596e-6*thrust_sp - 2.669e-8*y[3]*d_p_sp - 2.669e-8*p_sp*yp[3] + 2.669e-8*p_sp*d_p_sp + 2.669e-8*y[4]*yp[4] - 2.669e-8*y[4]*d_q_sp - 2.669e-8*q_sp*yp[4] + 2.669e-8*q_sp*d_q_sp + 5.339e-6*y[3]*y[6] - 5.339e-6*y[3]*int_psp - 5.339e-6*p_sp*y[6] + 5.339e-6*p_sp*int_psp + 5.339e-6*y[4]*y[7] - 5.339e-6*y[4]*int_qsp - 5.339e-6*q_sp*y[7] + 5.339e-6*q_sp*int_qsp + 3.423e-7*y[5]*y[8] - 3.423e-7*y[5]*int_rsp - 3.423e-7*r_sp*y[8] + 3.423e-7*r_sp*int_rsp + /*1.335e-10*yp[3]^2 + 1.335e-10*d_p_sp^2 + 1.335e-10*yp[4]^2 + 1.335e-10*d_q_sp^2 +*/ 5.339e-6*y[6]*y[6] + 5.339e-6*int_psp*int_psp + 5.339e-6*y[7]*y[7] + 5.339e-6*int_qsp*int_qsp + 2.382e-8*y[8]*y[8] + 2.382e-8*int_rsp*int_rsp - 2.669e-6*y[3]*p_sp - 2.669e-6*y[4]*q_sp - 2.46e-6*y[5]*r_sp + 1.335e-6*y[3]*y[3] + 1.335e-6*p_sp*p_sp + 1.335e-6*y[4]*y[4] + 1.335e-6*q_sp*q_sp + 1.23e-6*y[5]*y[5] + 1.23e-6*r_sp*r_sp + 8.542e-11*thrust_sp*thrust_sp /*- 2.669e-10*yp[3]*d_p_sp - 2.669e-10*yp[4]*d_q_sp*/ + 5.339e-8*yp[3]*y[6] - 5.339e-8*yp[3]*int_psp - 5.339e-8*d_p_sp*y[6] + 5.339e-8*d_p_sp*int_psp + 5.339e-8*yp[4]*y[7] - 5.339e-8*yp[4]*int_qsp - 5.339e-8*d_q_sp*y[7] + 5.339e-8*d_q_sp*int_qsp - 1.068e-5*y[6]*int_psp - 1.068e-5*y[7]*int_qsp - 4.764e-8*y[8]*int_rsp + 2.669e-8*y[3]*yp[3] + 0.007455;
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // Translation in Z
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              yp[13] = y[11] + 2*y[12] - 2*z_sp;
              
              //Time variable
              yp[14] = 1.0;
          } else if (syschoice == 13)
          {
              /* Z and roll pitch yaw rate control with a PI controller */
              static const double p_sp =0.0;
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const AAF Ct = interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_p = 0.0;
              static const AAF Ip_pq = 0.0;
              static const AAF Ip_pr = 0.0;
              static const AAF Ip_q = 0.0;
              static const AAF Ip_qr = interval(-1.04880447793, -1.03580464787);
              static const AAF Ip_r = 0.0;
              
              static const AAF Iq_p = 0.0;
              static const AAF Iq_pq = 0.0;
              static const AAF Iq_pr = interval(1.03470095927, 1.04749270535);
              static const AAF Iq_q = 0.0;
              static const AAF Iq_qr = 0.0;
              static const AAF Iq_r = 0.0;
              
              static const AAF Ir_p = 0.0;
              static const AAF Ir_pq = interval(-0.0162919189567, -0.0120891632629);
              static const AAF Ir_pr = 0.0;
              static const AAF Ir_q = 0.0;
              static const AAF Ir_qr = 0.0;
              static const AAF Ir_r = 0.0;
              
              static const AAF Im_xx = interval(71484.0524534, 71885.7226787);
              static const AAF Im_xy = 0.0;
              static const AAF Im_xz = 0.0;
              static const AAF Im_yy = interval(69441.6509547, 69834.7034512);
              static const AAF Im_yz = 0.0;
              static const AAF Im_zz = interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              // auto cosYaw = cos(y[2]*M_PI/180.0);
              // auto sinYaw = sin(y[2]*M_PI/180.0);
              
              auto tanPitch = tan(y[1]);
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto err_z = z_sp - y[12];
              
              auto velZ_sp = Kp_Z * err_z;
              auto thrust_Raw = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
              
              auto thrust = 1000.0*thrust_Raw + 36000;
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_p;
              yp[7] = err_q;
              yp[8] = err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // Z derivative coordinate
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              
              // Z integrale for thrust setpoint calculation
              yp[13] = velZ_sp - y[11];
          } else if(syschoice == 14)
          {
              /*Control in x , y , z yaw PI */
              
              // static const double roll_sp =10.0;
              // static const double pitch_sp = 0.0;
              // static const double yaw_sp = 0.0;
              static const double x_sp = 0.1;
              static const double y_sp = 0.0;
              static const double z_sp = 1.0;
              static const double yaw_sp = 0.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_X = 2.0;
              static const double Kp_VX = 25.0;
              static const double Ki_VX = 1.0;
              
              static const double Kp_Y = 2.0;
              static const double Kp_VY = 25.0;
              static const double Ki_VY = 1.0;
              
              static const double Kp_r = 6.0;
              static const double Kp_p = 6.0;
              static const double Kp_y = 6.0;
              
              static const double Ki_r = 3.0;
              static const double Ki_p = 3.0;
              static const double Ki_y = 1.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;//120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              // static const double Kp_Z = 2.0;
              // static const double Kp_VZ = 25.0;
              // static const double Ki_VZ = 15.0;
              
              // static const double Kp_X = 2.0;
              // static const double Kp_VX = 10.0;
              // static const double Ki_VX = 1.0;
              
              // static const double Kp_Y = 2.0;
              // static const double Kp_VY = 10.0;
              // static const double Ki_VY = 1.0;
              
              // static const double Kp_r = 6.0;
              // static const double Kp_p = 6.0;
              // static const double Kp_y = 6.0;
              
              // static const double Ki_r = 3.0;
              // static const double Ki_p = 3.0;
              // static const double Ki_y = 0.0;
              
              // static const double Kp_rr = 800.0;
              // static const double Kp_pr = 800.0;
              // static const double Kp_yr = 1000;//120.0;
              
              // static const double Ki_rr = 0.0;
              // static const double Ki_pr = 0.0;
              // static const double Ki_yr = 0.0;
              
              // static const double Kp_r = 3.5;
              // static const double Kp_p = 3.5;
              // static const double Kp_y = 6.0;
              
              // static const double Ki_r = 2.0;
              // static const double Ki_p = 2.0;
              // static const double Ki_y = 1.0;
              
              // static const double Kp_rr = 70.0;
              // static const double Kp_pr = 70.0;
              // static const double Kp_yr = 70.0;
              
              // static const double Ki_rr = 0.0;
              // static const double Ki_pr = 0.0;
              // static const double Ki_yr = 16.7;
              
              
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              // Error in each x,y,z target
              auto err_x = x_sp - y[12];
              auto err_y = y_sp - y[13];
              auto err_z = z_sp - y[14];
              
              // I remove the integrator part for simplicity (reduce number of variables)
              // In the source code the derivative coefficient are zero so the first PID is just a Proportionnal here
              auto velX_sp = Kp_X * err_x;
              auto velY_sp = Kp_Y * err_y;
              auto velZ_sp = Kp_Z * err_z;
              
              // Second PID is taking as target the velocity output from first PID
              // And it's siimply a PI here. There is actually a Kd coefficient in the
              // crazyflie source code but I assume it's 0 here due for some reason I mentionned
              // In the report
              auto roll_sp_raw    = Kp_VX * (velX_sp - y[9]) + Ki_VX * y[15];
              auto pitch_sp_raw   = Kp_VY * (velY_sp - y[10]) + Ki_VY * y[16];
              auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[17];
              
              // The value obtained above doesn't take in account the yaw value -> correction
              auto roll_sp =  -cosYaw * pitch_sp_raw + sinYaw * roll_sp_raw ;
              auto pitch_sp = cosYaw * roll_sp_raw + sinYaw * pitch_sp_raw ;
              auto thrust = 1000.0*thrust_Raw + 36000;
              
              // Here is the euler angle PID loop control
              auto err_roll =  roll_sp - y[0];
              auto err_pitch = pitch_sp - y[1];
              auto err_yaw = yaw_sp - y[2];
              
              auto p_sp = Kp_r * (err_roll) + Ki_r * y[6];//y[18];
              auto q_sp = Kp_p * (err_pitch) + Ki_p * y[7];//y[19];
              auto r_sp = Kp_y * (err_yaw) + Ki_y * y[8];//y[20];
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto cmd_r = y[6]*Ki_rr +  err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr +  err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr +  err_r*Kp_yr;
              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_roll;//err_p;
              yp[7] = err_pitch;//err_q;
              yp[8] = err_yaw;//err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              //std::cout << "u : " << getAAF(pitch_sp).convert_int() << std::endl;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // X , Y ,Z derivative coordinate
              yp[12] = y[11]*(sinRoll*sinYaw + cosRoll*cosYaw*sinPitch) - y[10]*(cosRoll*sinYaw - cosYaw*sinPitch*sinRoll) + cosPitch*cosYaw*y[9];;
              yp[13] = y[10]*(cosRoll*cosYaw + sinPitch*sinRoll*sinYaw) - y[11]*(cosYaw*sinRoll - cosRoll*sinPitch*sinYaw) + cosPitch*sinYaw*y[9];;
              yp[14] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // Z integrale for thrust setpoint calculation
              yp[15] = velX_sp - y[9];
              yp[16] = velY_sp - y[10];
              yp[17] = velZ_sp - y[11];
              
              //yp[18] = err_roll;
              //yp[19] = err_pitch;
              //yp[20] = err_yaw;
              //std::cout << getAAF(thrust).convert_int() << endl;
          } else if (syschoice == 15) // Mellinger controller
          {
              //static const double x_sp = 0.1;
              //static const double y_sp = 0.1;
              //static const double z_sp = 1.0;
              //static const double yaw_sp = 0.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double kp_xy = 0.4;
              static const double kd_xy = 0.2;
              static const double ki_xy = 0.05;
              
              static const double kp_z = 1.25;
              static const double kd_z = 0.4;
              static const double ki_z = 0.05;
              
              static const double kR_xy = 70000;
              static const double kw_xy = 20000;
              static const double ki_m_xy = 0.0;
              
              static const double kR_z = 60000;
              static const double kw_z = 12000;
              static const double ki_m_z = 500;
              
              static const double kd_omega_rp = 200;
              
              static const double m_cf = 0.032;
              static const double massThrust = 132000;
              static const double g = 9.81;
              
              // We maintain the yaw angle to 0
              //static const double cosYaw = 1.0;
              //static const double sinYaw = 0.0;
              
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;
              
              static const AAF Im_xx = 1.0/Ixx;
              static const AAF Im_yy = 1.0/Iyy;
              static const AAF Im_zz = 1.0/Izz;
              
              // setpoint values
              auto x_sp = 0.1*y[19];
              auto vx_sp = 0.1;
              auto ax_sp = 0.0;
              
              auto y_sp = 0.0;//0.1*y[19];
              auto vy_sp = 0.0;
              auto ay_sp = 0.0;
              
              auto z_sp = 1.0;
              auto vz_sp = 0.0;
              auto az_sp = 0.0;
              
              auto yaw_sp = 0.0; // In radian
              auto gyro_x = 0.0;
              auto gyro_y = 0.0;
              auto gyro_z = 0.0;
              
              auto s = 2.0/(y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3]);
              //auto s = 2.0;
              
              auto err_x = x_sp - y[7];
              auto err_vx = vx_sp - y[10];
              auto err_y = y_sp - y[8];
              auto err_vy = vy_sp - y[11];
              auto err_z = z_sp - y[9];
              auto err_vz = vz_sp - y[12];
              
              auto temp1 = s*(y[1]*y[1] + y[2]*y[2]) - 1;
              auto temp2 = s*(y[1]*y[1] + y[3]*y[3]) - 1;
              auto temp3 = s*(y[2]*y[2] + y[3]*y[3]) - 1;
              auto temp4 = s*(y[0]*y[2] + y[1]*y[3]);
              auto temp5 = s*(y[0]*y[2] - y[1]*y[3]);
              auto temp6 = s*(y[0]*y[1] - y[2]*y[3]);
              auto temp7 = s*(y[0]*y[1] + y[2]*y[3]);
              auto temp8 = s*(y[0]*y[3] - y[1]*y[2]);
              auto temp9 = s*(y[0]*y[3] + y[1]*y[2]);
              
              auto z_axis_x =  temp4;//s*(y[1]*y[3]+y[2]*y[0]) ;
              auto z_axis_y = -temp6;//s*(y[2]*y[3]-y[1]*y[0]) ;
              auto z_axis_z = -temp1;//1.0 - s*(y[1]*y[1] + y[2]*y[2]);
              
              auto target_thrust_x = m_cf * ax_sp + kp_xy*err_x + kd_xy*err_vx + ki_xy * y[13];
              auto target_thrust_y = m_cf * ay_sp + kp_xy*err_y + kd_xy*err_vy + ki_xy * y[14];
              auto target_thrust_z = m_cf * (az_sp + g) + kp_z*err_z + kd_z*err_vz + ki_z * y[15];
              
              auto thrust = massThrust*(target_thrust_x * z_axis_x + target_thrust_y*z_axis_y + target_thrust_z*z_axis_z);
              
           //   cout << getAAF(thrust).convert_int() << endl;
              auto target_norm_inv = 1.0/sqrt(target_thrust_x*target_thrust_x + target_thrust_y*target_thrust_y + target_thrust_z*target_thrust_z);
              auto z_axis_desired_x = target_thrust_x * target_norm_inv;
              auto z_axis_desired_y = target_thrust_y * target_norm_inv;
              auto z_axis_desired_z = target_thrust_z * target_norm_inv;
              
              // With the yaw setpoint set always to 0
              auto y_axis_desired_x = 0.0;
              auto y_axis_desired_y = z_axis_desired_z ;
              auto y_axis_desired_z = -z_axis_desired_y;
              auto y_axis_norm = 1.0/sqrt(y_axis_desired_y*y_axis_desired_y + y_axis_desired_z*y_axis_desired_z);
              y_axis_desired_y = y_axis_desired_y * y_axis_norm;
              y_axis_desired_z = y_axis_desired_z * y_axis_norm;
              
              auto x_axis_desired_x = y_axis_desired_y*z_axis_desired_z - y_axis_desired_z*z_axis_desired_y;
              auto x_axis_desired_y = y_axis_desired_z*z_axis_desired_x; //- y_axis_desired_x*z_axis_desired_z;
              auto x_axis_desired_z = /*y_axis_desired_x*z_axis_desired_y*/ - y_axis_desired_y*z_axis_desired_x;
              
              auto eR_x = y_axis_desired_z*(temp1) - z_axis_desired_y*(temp2) /*- y_axis_desired_x*temp4*/ + y_axis_desired_y*temp6 - z_axis_desired_x*temp8 + z_axis_desired_z*temp7;
              auto eR_y = z_axis_desired_x*(temp3) - x_axis_desired_z*(temp1) + x_axis_desired_x*temp4 - x_axis_desired_y*temp6 - z_axis_desired_y*temp9 + z_axis_desired_z*temp5;
              auto eR_z = x_axis_desired_y*(temp2) /*- y_axis_desired_x*(temp3)*/ + x_axis_desired_x*temp8 - x_axis_desired_z*temp7 + y_axis_desired_y*temp9 - y_axis_desired_z*temp5;
              
              auto cmd_r = -kR_xy * eR_x + kw_xy * (gyro_x - y[4]) + ki_m_xy * y[16];
              auto cmd_p = -kR_xy * eR_y + kw_xy * (gyro_y - y[5]) + ki_m_xy * y[17];
              auto cmd_y = -kR_z * eR_z + kw_z * (gyro_z - y[6]) + ki_m_z * y[18];
              
              //cout << getAAF(cmd_r).convert_int() << endl;
              //cout << getAAF(cmd_p).convert_int() << endl;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Quaternion dynamics
              yp[0] = -0.5 * (y[1]*y[4] + y[2]*y[5] + y[3]*y[6]);
              yp[1] =  0.5 * (y[0]*y[4] - y[3]*y[5] + y[2]*y[6]);
              yp[2] =  0.5 * (y[3]*y[4] + y[0]*y[5] - y[1]*y[6]);
              yp[3] = -0.5 * (y[2]*y[4] - y[1]*y[5] - y[0]*y[6]);
              
              // p q r rotation dynamics
              yp[4] = Ip_qr * y[5] * y[6] + Im_xx * Mx ;
              yp[5] = Iq_pr * y[4] * y[6] + Im_yy * My ;
              yp[6] = Ir_pq * y[4] * y[5] + Im_zz * Mz ;
              
              // Position x , y , z in world frame
              yp[7] = y[10];
              yp[8] = y[11];
              yp[9] = y[12];
              
              //Velocity vx vy and vz in world frame
              auto F_m =  F/m_cf;
              yp[10] = z_axis_x * F_m;
              yp[11] = z_axis_y * F_m;
              yp[12] = z_axis_z * F_m - g;
              
              // integrale of error in x y and z
              yp[13] = err_x;
              yp[14] = err_y;
              yp[15] = err_z;
              
              // integrale in error for orientation
              yp[16] = -eR_x;
              yp[17] = -eR_y;
              yp[18] = -eR_z;
              
              // Derivative of time
              yp[19] = 1.0;
          } else if (syschoice == 16) // Mellinger controller
          {
              /* Euler anglle not quaternion */
              //static const double x_sp = 0.1;
              //static const double y_sp = 0.1;
              //static const double z_sp = 1.0;
              //static const double yaw_sp = 0.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double kp_xy = 0.4;
              static const double kd_xy = 0.2;
              static const double ki_xy = 0.05;
              
              static const double kp_z = 1.25;
              static const double kd_z = 0.4;
              static const double ki_z = 0.05;
              
              static const double kR_xy = 70000;
              static const double kw_xy = 20000;
              static const double ki_m_xy = 0.0;
              
              static const double kR_z = 60000;
              static const double kw_z = 12000;
              static const double ki_m_z = 500;
              
              static const double kd_omega_rp = 200;
              
              static const double m_cf = 0.032;
              static const double massThrust = 132000;
              static const double g = 9.81;
              
              // We maintain the yaw angle to 0
              //static const double cosYaw = 1.0;
              //static const double sinYaw = 0.0;
              
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;
              
              static const AAF Im_xx = 1.0/Ixx;
              static const AAF Im_yy = 1.0/Iyy;
              static const AAF Im_zz = 1.0/Izz;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              //yp[15] = y[20]*(sinRoll*sinYaw + cosRoll*cosYaw*sinPitch) - y[19]*(cosRoll*sinYaw - cosYaw*sinPitch*sinRoll) + cosPitch*cosYaw*y[18];
              //yp[16] = y[19]*(cosRoll*cosYaw + sinPitch*sinRoll*sinYaw) - y[20]*(cosYaw*sinRoll - cosRoll*sinPitch*sinYaw) + cosPitch*sinYaw*y[18];
              //yp[17] = cosPitch*cosRoll*y[20] - sinPitch*y[18] + cosPitch*sinRoll*y[19];
              
              auto temp1 = cosPitch*cosYaw;
              auto temp2 = (-cosRoll*sinYaw + cosYaw*sinPitch*sinRoll);
              auto temp3 = (sinRoll*sinYaw + cosRoll*cosYaw*sinPitch);
              auto temp4 = cosPitch*sinYaw;
              auto temp5 = (cosRoll*cosYaw + sinPitch*sinRoll*sinYaw);
              auto temp6 = (-cosYaw*sinRoll + cosRoll*sinPitch*sinYaw);
              auto temp7 = - sinPitch;
              auto temp8 = cosPitch*sinRoll;
              auto temp9 = cosPitch*cosRoll;
              
              // setpoint values
              auto x_sp = 0.1*y[18];
              auto vx_sp = 0.1;
              auto ax_sp = 0.0;
              
              auto y_sp = 0.1*y[18];
              auto vy_sp = 0.1;
              auto ay_sp = 0.0;
              
              auto z_sp = 1.0;
              auto vz_sp = 0.0;
              auto az_sp = 0.0;
              
              auto yaw_sp = 0.0;
              auto gyro_x = 0.0;
              auto gyro_y = 0.0;
              auto gyro_z = 0.0;
              
              auto err_x = x_sp - y[6];
              auto err_vx = vx_sp - y[9];
              auto err_y = y_sp - y[7];
              auto err_vy = vy_sp - y[10];
              auto err_z = z_sp - y[8];
              auto err_vz = vz_sp - y[11];
              
              
              auto target_thrust_x = m_cf * ax_sp + kp_xy*err_x + kd_xy*err_vx + ki_xy * y[12];
              auto target_thrust_y = m_cf * ay_sp + kp_xy*err_y + kd_xy*err_vy + ki_xy * y[13];
              auto target_thrust_z = m_cf * (az_sp + g) + kp_z*err_z + kd_z*err_vz + ki_z * y[14];
              
              auto thrust = massThrust*(target_thrust_x * temp7 + target_thrust_y*temp8 + target_thrust_z*temp9);
              
              auto inv_Fdez = 1.0/target_thrust_z;
              
              auto pitch_desired = atan(target_thrust_x * inv_Fdez);
              auto cosPitchDesired = cos(pitch_desired);
              auto sinPitchDesired = sin(pitch_desired);
              
              auto rollDesired = atan(-cosPitchDesired * target_thrust_y * inv_Fdez);
              auto cosRollDesired = cos(rollDesired);
              auto sinRollDesired = sin(rollDesired);
              
              auto x_axis_desired_x = cosPitchDesired;
              auto x_axis_desired_y = 0.0;
              auto x_axis_desired_z = -sinPitchDesired;
              
              auto y_axis_desired_x = sinPitchDesired*sinRollDesired;
              auto y_axis_desired_y = cosRollDesired;
              auto y_axis_desired_z = sinRollDesired*cosPitchDesired;
              
              auto z_axis_desired_x = cosRollDesired*sinPitchDesired;
              auto z_axis_desired_y = -sinRollDesired;
              auto z_axis_desired_z = cosRollDesired*cosPitchDesired;
              
              auto eR_x = y_axis_desired_z*(-temp9) + z_axis_desired_y*(temp5) - y_axis_desired_x*temp3 - y_axis_desired_y*temp6 + z_axis_desired_x*temp2 + z_axis_desired_z*temp8;
              auto eR_y = z_axis_desired_x*(-temp1) + x_axis_desired_z*(temp9) + x_axis_desired_x*temp3 + x_axis_desired_y*temp6 - z_axis_desired_y*temp4 - z_axis_desired_z*temp7;
              auto eR_z = x_axis_desired_y*(-temp5) + y_axis_desired_x*(temp1) - x_axis_desired_x*temp2 - x_axis_desired_z*temp8 + y_axis_desired_y*temp4 + y_axis_desired_z*temp7;
              
              auto cmd_r = -kR_xy * eR_x + kw_xy * (gyro_x - y[3]) + ki_m_xy * y[15];
              auto cmd_p = -kR_xy * eR_y + kw_xy * (gyro_y - y[4]) + ki_m_xy * y[16];
              auto cmd_y = -kR_z * eR_z + kw_z * (gyro_z - y[5]) + ki_m_z * y[17];
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Quaternion dynamics
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p q r rotation dynamics
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // Position x , y , z in world frame
              yp[6] = y[11]*(sinRoll*sinYaw + cosRoll*cosYaw*sinPitch) - y[10]*(cosRoll*sinYaw - cosYaw*sinPitch*sinRoll) + cosPitch*cosYaw*y[9];;
              yp[7] = y[10]*(cosRoll*cosYaw + sinPitch*sinRoll*sinYaw) - y[11]*(cosYaw*sinRoll - cosRoll*sinPitch*sinYaw) + cosPitch*sinYaw*y[9];;
              yp[8] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              
              //Velocity vx vy and vz in world frame
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m_cf - g*cosPitch*cosRoll;
              
              // integrale of error in x y and z
              yp[12] = err_x;
              yp[13] = err_y;
              yp[14] = err_z;
              
              // integrale in error for orientation
              yp[15] = -eR_x;
              yp[16] = -eR_y;
              yp[17] = -eR_z;
              
              // Derivative of time
              yp[18] = 1.0;
          } else if (syschoice == 17) // Mellinger controller
          {
              //static const double x_sp = 0.1;
              //static const double y_sp = 0.1;
              //static const double z_sp = 1.0;
              //static const double yaw_sp = 0.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double kp_xy = 0.4;
              static const double kd_xy = 0.2;
              static const double ki_xy = 0.05;
              
              static const double kp_z = 1.25;
              static const double kd_z = 0.4;
              static const double ki_z = 0.05;
              
              //static const double kR_xy = 70000;
              static const double kw_xy = 20000;
              //static const double ki_m_xy = 0.0;
              
              //static const double kR_z = 60000;
              static const double kw_z = 12000;
              //static const double ki_m_z = 500;
              
              
              static const double Kp_y = 6.0;
              
              //static const double kd_omega_rp = 200;
              
              static const double m_cf = 0.032;
              static const double massThrust = 132000;
              static const double g = 9.81;
              
              // We maintain the yaw angle to 0
              //static const double cosYaw = 1.0;
              //static const double sinYaw = 0.0;
              
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;
              
              static const AAF Im_xx = 1.0/Ixx;
              static const AAF Im_yy = 1.0/Iyy;
              static const AAF Im_zz = 1.0/Izz;
              
              // setpoint values
              auto x_sp = 0.1*y[15];
              auto vx_sp = 0.1;
              auto ax_sp = 0.0;
              
              auto y_sp = 0.0;//0.1*y[15];
              auto vy_sp = 0.0;
              auto ay_sp = 0.0;
              
              auto z_sp = 1.0;
              auto vz_sp = 0.0;
              auto az_sp = 0.0;
              
              auto yaw_sp = 0.0;
              auto gyro_x = 0.0;
              auto gyro_y = 0.0;
              auto gyro_z = 0.0;
              
              auto err_x = x_sp - y[6];
              auto err_vx = vx_sp - y[9];
              auto err_y = y_sp - y[7];
              auto err_vy = vy_sp - y[10];
              auto err_z = z_sp - y[8];
              auto err_vz = vz_sp - y[11];
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              auto temp1 = cosPitch*cosYaw;
              auto temp2 = (-cosRoll*sinYaw + cosYaw*sinPitch*sinRoll);
              auto temp3 = (sinRoll*sinYaw + cosRoll*cosYaw*sinPitch);
              auto temp4 = cosPitch*sinYaw;
              auto temp5 = (cosRoll*cosYaw + sinPitch*sinRoll*sinYaw);
              auto temp6 = (-cosYaw*sinRoll + cosRoll*sinPitch*sinYaw);
              auto temp7 = - sinPitch;
              auto temp8 = cosPitch*sinRoll;
              auto temp9 = cosPitch*cosRoll;
              
              auto target_thrust_x = m_cf * ax_sp + kp_xy*err_x + kd_xy*err_vx + ki_xy * y[12];
              auto target_thrust_y = m_cf * ay_sp + kp_xy*err_y + kd_xy*err_vy + ki_xy * y[13];
              auto target_thrust_z = m_cf * (az_sp + g) + kp_z*err_z + kd_z*err_vz + ki_z * y[14];
              
              auto thrust = massThrust * (target_thrust_x * temp3 + target_thrust_y*temp6 + target_thrust_z*temp9);
              
              // auto Fzy_sqr = target_thrust_z*target_thrust_z + target_thrust_y*target_thrust_y;
              
              // auto target_norm_inv = 1.0/sqrt(Fzy_sqr + target_thrust_x*target_thrust_x);
              // auto z_axis_desired_x = target_thrust_x * target_norm_inv;
              // auto z_axis_desired_y = target_thrust_y * target_norm_inv;
              // auto z_axis_desired_z = target_thrust_z * target_norm_inv;
              
              // // With the yaw setpoint set always to 0
              // //auto y_axis_desired_x = 0.0;
              // auto y_axis_norm = 1.0/sqrt(Fzy_sqr);
              // auto y_axis_desired_y = target_thrust_z*y_axis_norm ;
              // auto y_axis_desired_z = -target_thrust_y*y_axis_norm;
              
              // auto x_axis_desired_x = y_axis_desired_y*z_axis_desired_z - y_axis_desired_z*z_axis_desired_y;
              // auto x_axis_desired_y = y_axis_desired_z*z_axis_desired_x; //- y_axis_desired_x*z_axis_desired_z;
              // auto x_axis_desired_z = /*y_axis_desired_x*z_axis_desired_y*/ - y_axis_desired_y*z_axis_desired_x;
              
              // auto eR_x = temp2*z_axis_desired_x - temp6*y_axis_desired_y - temp9*y_axis_desired_z /*- temp3*y_axis_desired_x*/ + temp5*z_axis_desired_y + temp8*z_axis_desired_z;
              // auto eR_y = temp3*x_axis_desired_x + temp6*x_axis_desired_y + temp9*x_axis_desired_z - temp1*z_axis_desired_x - temp4*z_axis_desired_y - temp7*z_axis_desired_z;
              // auto eR_z = /*temp1*y_axis_desired_x*/ - temp5*x_axis_desired_y - temp8*x_axis_desired_z - temp2*x_axis_desired_x + temp4*y_axis_desired_y + temp7*y_axis_desired_z;*/
              
              //auto eR_x = -(target_thrust_y + sinRoll*target_norm_inv);
              //auto eR_y = target_thrust_x - cosRoll*sinPitch*target_norm_inv;
              //auto eR_z = target_thrust_z*target_norm_inv - temp9;
              
              //auto target_wz = Kp_y * ()
              //cout << getAAF(eR_x).convert_int() << endl;
              //cout << getAAF(eR_y).convert_int() << endl;
              //cout << getAAF(eR_z).convert_int() << endl;
              
              auto target_wx =  -(target_thrust_x * temp2 + target_thrust_y * temp5 + target_thrust_z * temp8);
              auto target_wy =   (target_thrust_x * temp1 + target_thrust_y * temp4 + target_thrust_z * temp7);
              
              auto cmd_r = /*-kR_xy * eR_x +*/ 2*kw_xy * (target_wx + gyro_x - y[3]) ;
              auto cmd_p = /*-kR_xy * eR_y +*/ 2*kw_xy * (target_wy + gyro_y - y[4]) ;
              auto cmd_y = /*-kR_z  * eR_z +*/ kw_z * (-6.0*y[2] + gyro_z - y[5]) ;
              //cout << getAAF(cmd_r).convert_int() << endl;
              //cout << getAAF(cmd_p).convert_int() << endl;
              
              //thrust *= massThrust;
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // euler angle dynamics
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              /*cout << getAAF(cosRoll).convert_int() << endl;
               cout << getAAF(sinRoll).convert_int() << endl;
               cout << getAAF(cosPitch) << endl;
               cout << getAAF(sinPitch) << endl;
               cout << "---------------" << endl;*/
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              //Position in world frame
              // yp[6] = y[11]*temp3 + y[10]*temp2 + y[9]*temp1;
              // yp[7] = y[11]*temp6 + y[10]*temp5 + y[9]*temp4;
              // yp[8] = y[11]*temp9 + y[10]*temp8 + y[9]*temp7;
              yp[6] = y[9];
              yp[7] = y[10];
              yp[8] = y[11];
              
              //Velocity u v and w in body frame
              // yp[9] =  y[5]*y[10] - y[4]*y[11] - g*temp7;
              // yp[10] = y[3]*y[11] - y[5]*y[9]  - g*temp8;
              // yp[11] = y[4]*y[9]  - y[3]*y[10] + F/m_cf - g*temp9;
              auto F_m = F/m_cf;
              yp[9] = temp3 * F_m;
              yp[10] = temp6 * F_m;
              yp[11] = temp9 * F_m - g;
              
              // integrale of error in x y and z
              yp[12] = err_x;
              yp[13] = err_y;
              yp[14] = err_z;
              
              // Derivative of time
              yp[15] = 1.0;
          }else if(syschoice == 18) // paper example
          {
              static const double p_sp = 1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              auto err_z = z_sp - y[12];
              
              auto velZ_sp = Kp_Z * err_z;
              
              //Z derivative coordinate
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // Z integrale for thrust setpoint calculation
              yp[13] = velZ_sp - yp[12];
              
              //auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
              auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
              
              auto thrust = 1000.0*thrust_Raw + 36000;
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_p;//err_p;
              yp[7] = err_q;//err_q;
              yp[8] = err_r;//err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // //Z derivative coordinate
              // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // // Z integrale for thrust setpoint calculation
              // yp[13] = velZ_sp - y[11];
          }
          else if (syschoice == 19) // crazyflie - paper example with  uncertain parameters
          {
              static const double p_sp = 1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              // constant parameters Ct and Cd
              yp[14] = 0; // Ct interval(1.28e-8 , 1.29e-8);
              yp[15] = 0; // Cd interval(7.64e-11 , 7.65e-11);
             // static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
             // static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              auto err_z = z_sp - y[12];
              
              auto velZ_sp = Kp_Z * err_z;
              
              //Z derivative coordinate
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // Z integrale for thrust setpoint calculation
              yp[13] = velZ_sp - yp[12];
              
              //auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
              auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
              
              auto thrust = 1000.0*thrust_Raw + 36000;
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              auto Mx = ((4*y[14]*d*thrust*C1*C1 + 4*C2*y[14]*d*C1)*cmd_r + (-4*C1*C1*y[14]*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*y[14]*d*cmd_r*cmd_y + (4*y[14]*d*thrust*C1*C1 + 4*C2*y[14]*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*y[15]*cmd_r*cmd_p + (8*y[15]*thrust*C1*C1 + 8*C2*y[15]*C1)*cmd_y);
              auto F  = y[14]*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*y[14]*C1*C2*thrust + 4*y[14]*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_p;//err_p;
              yp[7] = err_q;//err_q;
              yp[8] = err_r;//err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // //Z derivative coordinate
              // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // // Z integrale for thrust setpoint calculation
              // yp[13] = velZ_sp - y[11];
          }
          else if (syschoice == 20) // crazyflie - paper example with more uncertain parameters
          {
              static const double p_sp = 1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
             
              
              // constant parameters Ct and Cd
              yp[14] = 0; // Ct interval(1.28e-8 , 1.29e-8);
              yp[15] = 0; // Cd interval(7.64e-11 , 7.65e-11);
              // static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              // static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              /* Crazyflie trajectory tracking article*/
              
              
              // constant parameters Ixx, Iyy, Izz
              yp[16] = 0;
              yp[17] = 0;
              yp[18] = 0;
              
              // valeurs ponctuelles et intervalles calculees par inversion pas coherentes, bizarre ?
           //   static const double Ixx = 1.657171e-5;
           //   static const double Iyy = 1.6655602e-5;
           //   static const double Izz = 2.9261652e-5;
              
              // Ixx = interval(0.00001391096817,0.00001398913416);
              // Iyy = interval(0.00001431952812,0.00001440057928);
              // Izz = interval(0.00002880845915,0.00002899182825);
              
              
              auto Im_xx = 1.0/y[16];//interval(71484.0524534, 71885.7226787);
              auto Im_yy = 1.0/y[17];//interval(69441.6509547, 69834.7034512);
              auto Im_zz = 1.0/y[18];//interval(34492.4780616, 34712.0265858);
              
              auto Ir_pq = (y[16]-y[17])/y[18];//interval(-0.0162919189567, -0.0120891632629);
              auto Ip_qr = (y[17]-y[18])/y[16];//interval(-1.04880447793, -1.03580464787);
              auto Iq_pr = (y[18]-y[16])/y[17];//interval(1.03470095927, 1.04749270535);
              
            
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              auto err_z = z_sp - y[12];
              
              auto velZ_sp = Kp_Z * err_z;
              
              //Z derivative coordinate
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // Z integrale for thrust setpoint calculation
              yp[13] = velZ_sp - yp[12];
              
              //auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
              auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
              
              auto thrust = 1000.0*thrust_Raw + 36000;
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              auto Mx = ((4*y[14]*d*thrust*C1*C1 + 4*C2*y[14]*d*C1)*cmd_r + (-4*C1*C1*y[14]*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*y[14]*d*cmd_r*cmd_y + (4*y[14]*d*thrust*C1*C1 + 4*C2*y[14]*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*y[15]*cmd_r*cmd_p + (8*y[15]*thrust*C1*C1 + 8*C2*y[15]*C1)*cmd_y);
              auto F  = y[14]*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*y[14]*C1*C2*thrust + 4*y[14]*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_p;//err_p;
              yp[7] = err_q;//err_q;
              yp[8] = err_r;//err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // //Z derivative coordinate
              // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // // Z integrale for thrust setpoint calculation
              // yp[13] = velZ_sp - y[11];
          }
          else if (syschoice == 21) // paper example but with additional x and y
          {
              static const double p_sp = 1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              
              static const double Ixx = 1.657171e-5;
              static const double Iyy = 1.6655602e-5;
              static const double Izz = 2.9261652e-5;
              
              static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
              static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
              
              static const double Kp_Z = 2.0;
              static const double Kp_VZ = 25.0;
              static const double Ki_VZ = 15.0;
              
              static const double Kp_rr = 250.0;
              static const double Kp_pr = 250.0;
              static const double Kp_yr = 120.0;
              
              static const double Ki_rr = 500.0;
              static const double Ki_pr = 500.0;
              static const double Ki_yr = 16.7;
              
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              auto err_z = z_sp - y[12];
              
              auto velZ_sp = Kp_Z * err_z;
              
              //Z derivative coordinate
              yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // Z integrale for thrust setpoint calculation
              yp[13] = velZ_sp - yp[12];
              
              // derivative x
              yp[14] = cosPitch*cosYaw*y[9]+(cosYaw*sinPitch*sinRoll-cosRoll*sinYaw)*y[10]+(sinRoll*sinYaw+cosRoll*cosYaw*sinPitch)*y[11];
              // derivative y
              yp[15] = cosPitch*sinYaw*y[9]+(cosRoll*cosYaw+sinPitch*sinRoll*sinYaw)*y[10]+(cosRoll*sinPitch*sinYaw-cosYaw*sinRoll)*y[11];
              
              //auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
              auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
              
              auto thrust = 1000.0*thrust_Raw + 36000;
              
              auto err_p = p_sp - y[3];
              auto err_q = q_sp - y[4];
              auto err_r = r_sp - y[5];
              
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
              
              // Roll , pitch , yaw derivatives
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              // integrale of error in p , q and r
              yp[6] = err_p;//err_p;
              yp[7] = err_q;//err_q;
              yp[8] = err_r;//err_r;
              
              // derivatives of body speed u , v and w
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // //Z derivative coordinate
              // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // // Z integrale for thrust setpoint calculation
              // yp[13] = velZ_sp - y[11];
          }
          else if (syschoice == 22) // Franck's system 19 (document p 12?) Roll pitch yaw z control example
              {
                  static const double roll_sp = 1.0 * M_PI / 180.0;
                  static const double pitch_sp = -1.0 * M_PI / 180.0;
                  static const double yaw_sp = 0.0;
                  
                  static const double z_sp = 1.0;
                  
                  static const double C1 = 0.04076521;
                  static const double C2 = 380.8359;
                  static const double d = 0.046/sqrt(2.0);
                  
                  static const double Ixx = 1.657171e-5;
                  static const double Iyy = 1.6655602e-5;
                  static const double Izz = 2.9261652e-5;
                  
                  static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
                  static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
                  
                  static const double Kp_Z = 2.0;
                  static const double Kp_VZ = 25.0;
                  static const double Ki_VZ = 15.0;
                  
                  static const double Kp_rr = 250.0;
                  static const double Kp_pr = 250.0;
                  static const double Kp_yr = 120.0;
                  
                  static const double Ki_rr = 500.0;
                  static const double Ki_pr = 500.0;
                  static const double Ki_yr = 16.7;
                  
                  static const double Kp_r = 6.0;
                  static const double Kp_p = 6.0;
                  static const double Kp_y = 6.0;
                  
                  static const double Ki_r = 3.0;
                  static const double Ki_p = 3.0;
                  static const double Ki_y = 1.0;
                  
                  /* Crazyflie trajectory tracking article*/
                  static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
                  static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
                  static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
                  
                  static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
                  static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
                  static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
                  
                  static const double g = 9.8;
                  static const double m = 0.028;
                  
                  auto cosRoll = cos(y[0]);
                  auto sinRoll = sin(y[0]);
                  
                  auto cosPitch = cos(y[1]);
                  auto sinPitch = sin(y[1]);
                  
                  auto cosYaw = cos(y[2]);
                  auto sinYaw = sin(y[2]);
                  
                  auto tanPitch = tan(y[1]);
                  
                  // integrale of error in roll, pitch , yaw
                  yp[14] = roll_sp - y[0];
                  yp[15] = pitch_sp - y[1];
                  yp[16] = yaw_sp - y[2];
                  
                  // From error above, p_sp, q_sp and r_sp are found
                  auto p_sp = Kp_r * yp[14] + Ki_r * y[14];
                  auto q_sp = Kp_p * yp[15] + Ki_p * y[15];;
                  auto r_sp = Kp_y * yp[16] + Ki_y * y[16];;
                  
                  auto err_z = z_sp - y[12];
                  auto velZ_sp = Kp_Z * err_z;
                  
                  //Z derivative coordinate
                  yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
                  // Z integrale for thrust setpoint calculation
                  yp[13] = velZ_sp - yp[12];
                  
                  //auto thrust_Raw     = Kp_VZ * (velZ_sp - y[11]) + Ki_VZ * y[13];
                  auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
                  
                  auto thrust = 1000.0*thrust_Raw + 36000;
                  
                  auto err_p = p_sp - y[3];
                  auto err_q = q_sp - y[4];
                  auto err_r = r_sp - y[5];
                  
                  auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr;
                  auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr;
                  auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr;
                  //std:cout << getAAF(cmd_p).convert_int() << std::endl;
                  
                  auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y);
                  auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p);
                  auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y);
                  auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2;
                  
                  // Roll , pitch , yaw derivatives
                  yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
                  yp[1] = y[4]*cosRoll - y[5]*sinRoll;
                  yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
                  
                  // p , q and r derivatives
                  yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
                  yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
                  yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
                  
                  // integrale of error in p , q and r
                  yp[6] = err_p;//err_p;
                  yp[7] = err_q;//err_q;
                  yp[8] = err_r;//err_r;
                  
                  // derivatives of body speed u , v and w
                  yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch;
                  yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll;
                  yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
                  
                  // //Z derivative coordinate
                  // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
                  // // Z integrale for thrust setpoint calculation
                  // yp[13] = velZ_sp - y[11];
              }
          else if (syschoice == 23) // benchmark quadcopter ARCH 2018
          {
      //        static const double p_sp = 0.0; //1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              
              static const double z_sp = 1.0;
              static const double roll_sp = 0.0; // desired roll
              static const double pitch_sp = 0.0; // desired pitch
              
              static const double M = 1;
              static const double R2 = 0.01;
              static const double l2 = 0.25;
              static const double Mrotor = 0.1;
              static const double m = M + 4*Mrotor;
              
              static const double Ixx = 2/5.*M*R2 + 2*l2*Mrotor;// 1.657171e-5;
              static const double Iyy = Ixx; //1.6655602e-5;
              static const double Izz = Ixx + 2*l2*Mrotor;// 2.9261652e-5;
              
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.81;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              // Z
              yp[9] = - cosPitch*cosRoll*y[8] + sinPitch*y[6] - cosPitch*sinRoll*y[7];
              //  x
              yp[10] = cosPitch*cosYaw*y[6]+(cosYaw*sinPitch*sinRoll-cosRoll*sinYaw)*y[7]+(sinRoll*sinYaw+cosRoll*cosYaw*sinPitch)*y[8];
              //  y
              yp[11] = cosPitch*sinYaw*y[6]+(cosRoll*cosYaw+sinPitch*sinRoll*sinYaw)*y[7]+(cosRoll*sinPitch*sinYaw-cosYaw*sinRoll)*y[8];
              
              auto Mx = -(y[0]-roll_sp) - y[3];
              auto My = -(y[1]-pitch_sp) - y[4];
              auto Mz = 0; // heading uncontrolled
              auto F  = m*g - 10*(y[9]-z_sp) + 3*y[8];
              
              // Roll , pitch , yaw
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives roll (rate pitch rate yaw rate)
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              //  body speed u , v and w (linear velocities)
              yp[6] = y[5]*y[7] - y[4]*y[8] - g*sinPitch;
              yp[7] = y[3]*y[8] - y[5]*y[6] + g*cosPitch*sinRoll;
              yp[8] = y[4]*y[6] - y[3]*y[7] - F/m + g*cosPitch*cosRoll;
              
          }
          else if (syschoice == 24) // benchmark quadcopter ARCH 2018 modified for initial conditions/setpoints of HSCC 2019
          {
          //    static const double p_sp = 0.0; //1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              
              static const double z_sp = 1.0;
              static const double roll_sp = 0.0; // desired roll
              static const double pitch_sp = 1.0*M_PI/180.0; // desired pitch
              
              static const double M = 1;
              static const double R2 = 0.01;
              static const double l2 = 0.25;
              static const double Mrotor = 0.1;
              static const double m = M + 4*Mrotor;
              
              static const double Ixx = 2/5.*M*R2 + 2*l2*Mrotor;// 1.657171e-5;
              static const double Iyy = Ixx; //1.6655602e-5;
              static const double Izz = Ixx + 2*l2*Mrotor;// 2.9261652e-5;
              
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.81;
              
              auto cosRoll = cos(y[0]);
              auto sinRoll = sin(y[0]);
              
              auto cosPitch = cos(y[1]);
              auto sinPitch = sin(y[1]);
              
              auto cosYaw = cos(y[2]);
              auto sinYaw = sin(y[2]);
              
              auto tanPitch = tan(y[1]);
              
              // Z
              yp[9] = - cosPitch*cosRoll*y[8] + sinPitch*y[6] - cosPitch*sinRoll*y[7];
              //  x
              yp[10] = cosPitch*cosYaw*y[6]+(cosYaw*sinPitch*sinRoll-cosRoll*sinYaw)*y[7]+(sinRoll*sinYaw+cosRoll*cosYaw*sinPitch)*y[8];
              //  y
              yp[11] = cosPitch*sinYaw*y[6]+(cosRoll*cosYaw+sinPitch*sinRoll*sinYaw)*y[7]+(cosRoll*sinPitch*sinYaw-cosYaw*sinRoll)*y[8];
              
              auto Mx = -(y[0]-roll_sp) - y[3];
              auto My = -(y[1]-pitch_sp) - y[4];
              auto Mz = 0; // heading uncontrolled
              auto F  = m*g - 10*(y[9]-z_sp) + 3*y[8];
              
              // Roll , pitch , yaw
              yp[0] = y[3] + (y[5]*cosRoll + y[4]*sinRoll)*tanPitch;
              yp[1] = y[4]*cosRoll - y[5]*sinRoll;
              yp[2] = (y[5]*cosRoll + y[4]*sinRoll)/cosPitch;
              
              // p , q and r derivatives roll (rate pitch rate yaw rate)
              yp[3] = Ip_qr * y[4] * y[5] + Im_xx * Mx ;
              yp[4] = Iq_pr * y[3] * y[5] + Im_yy * My ;
              yp[5] = Ir_pq * y[3] * y[4] + Im_zz * Mz ;
              
              //  body speed u , v and w (linear velocities)
              yp[6] = y[5]*y[7] - y[4]*y[8] - g*sinPitch;
              yp[7] = y[3]*y[8] - y[5]*y[6] + g*cosPitch*sinRoll;
              yp[8] = y[4]*y[6] - y[3]*y[7] - F/m + g*cosPitch*cosRoll;
              
          }
          else if (syschoice == 25) // damaged aircraft (Control Oriented Learning on the Fly) - control constant
          {
              yp[0] = -0.021*y[0] + 0.122*y[1] - 0.322*y[2] + 0.01*y[5] + y[6];
              yp[1] = -0.209*y[0] - 0.53*y[1] = 2.21*y[2] - 0.064*y[5] - 0.044*y[6];
              yp[2] = 0.017*y[0] - 0.164*y[1] - 0.421*y[2] - 0.378*y[5] + 0.544*y[6] + 0.01*cos(y[0])*y[0] + 0.15*sin(y[0])*y[1]+0.5*sin(y[1])*y[6];
              yp[3] = y[2];
              yp[4] = - y[1] + 2.21*y[3];
              yp[5] = 0; // control u1, constant for now
              yp[6] = 0; // control u2, constant for now
          }
          else if (syschoice == 26) // damaged aircraft (Control Oriented Learning on the Fly) - control variable
          {
              yp[0] = -0.021*y[0] + 0.122*y[1] - 0.322*y[2] + 0.01*y[5] + y[6];
              yp[1] = -0.209*y[0] - 0.53*y[1] = 2.21*y[2] - 0.064*y[5] - 0.044*y[6];
              yp[2] = 0.017*y[0] - 0.164*y[1] - 0.421*y[2] - 0.378*y[5] + 0.544*y[6] + 0.01*cos(y[0])*y[0] + 0.15*sin(y[0])*y[1]+0.5*sin(y[1])*y[6];
              yp[3] = y[2];
              yp[4] = - y[1] + 2.21*y[3];
              yp[5] =  interval(-30.,30.); //AAF(interval(-2.,2.)) + y[4]; ; // control u1, constant for now
              yp[6] =  interval(-1.,1.); //AAF(interval(-2.,2.)) + y[4]; ; // control u2, constant for now
          }
    }
};


// define here  your DEE system yp = \dot y = f(y, y_prev, params)
class DdeFunc {
public:
    template <class C>
    void operator()(vector<C> &yp, vector<C> &y, vector<C> &y_prev, vector<AAF> &param_inputs) {
        if (syschoice == 1)  // running example
        {
           yp[0] = - y[0] * y_prev[0];
        }
        else if (syschoice == 2) //  example 5.15
        {
            yp[0] = y_prev[0] + y[1];
            yp[1] = y[0] - y_prev[0];
        }
        else if (syschoice == 3) // Xue et al. 2017 ex 3
        {
            yp[0] = 1.4*y[2] - 0.9*y_prev[0];
            yp[1] = 2.5*y[4] - 1.5*y[1];
            yp[2] = 0.6*y[6] - 0.8*y[2]*y[1];
            yp[3] = 2.0 - 1.3*y[3]*y[2];
            yp[4] = 0.7*y[0] - y[3]*y[4];
            yp[5] = 0.3*y[0] - 3.1*y[5];
            yp[6] = 1.8*y[5] - 1.5*y[6]*y[1];
        }
        else if (syschoice == 4) // Szczelina_1 2014
        {
            yp[0] = -y[0] + 0.5*y_prev[0]*y_prev[0];
        }
        else if (syschoice == 5) // Szczelina_2 2014
        {
            yp[0] = -y[0] -3.2 * y_prev[0] + y_prev[0]*y_prev[0]*y_prev[0];
        }
        else if (syschoice == 6)  // self-driving car
        {
            yp[0] = y[1];
            yp[1] = -params[0] *(y_prev[0] - 1.0) - params[1]*y_prev[1];   // pr = 1 is the reference position
        }
        else if (syschoice == 7)  // self-driving car with uncertain (but constant) coefficients
        {
            yp[0] = y[1];
            yp[1] = -y[2] *(y_prev[0] - 1.0) - y[3]*y_prev[1];   // pr = 1 is the reference position
            yp[2] = 0;
            yp[3] = 0;
        }
        else if (syschoice == 8)  // self-driving car with uncertain (but constant) coefficients
        {
            yp[0] = y[1];
            yp[1] = -param_inputs[2] *(y_prev[0] - 1.0) - param_inputs[3]*y_prev[1];   // pr = 1 is the reference position
        }
        else if (syschoice == 9) // Ex 4 of Zou CAV 2015
        {
            yp[0] = -y_prev[0]*y_prev[0]*y_prev[0];
        }
        else if (syschoice == 10) // platoon 5 vehicles
        {
            double a = 2.5;
            yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
            yp[1] = y[2]; //  'x2' : 'dx2',
            yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.-y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
            yp[3] = y[4]; //  'x3' : 'dx3',
            yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
            yp[5] = y[6]; //   'x4' : 'dx4',
            yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
            yp[7] = y[8];//  'x5' : 'dx5',
            yp[8] = a*(y_prev[6]-y[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
        }
        else if (syschoice == 11) // platoon 10 vehicles
        {
            double a = 2.5;
            yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
            yp[1] = y[2]; //  'x2' : 'dx2',
            yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.-y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
            yp[3] = y[4]; //  'x3' : 'dx3',
            yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
            yp[5] = y[6]; //   'x4' : 'dx4',
            yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
            yp[7] = y[8];//  'x5' : 'dx5',
            yp[8] = a*(y_prev[6]-y_prev[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
            yp[9] = y[10]; //   'x6' : 'dx7',
            yp[10] = a*(y_prev[8]-y_prev[10]);//   'dx6'
            yp[11] = y[12];//  'x7' : 'dx7',
            yp[12] = a*(y_prev[10]-y_prev[12]);//  'dx7' :
            yp[13] = y[14]; //   'x6' : 'dx7',
            yp[14] = a*(y_prev[12]-y_prev[14]);//   'dx6'
            yp[15] = y[16];//  'x7' : 'dx7',
            yp[16] = a*(y_prev[14]-y[16]);//  'dx7' :
            yp[17] = y[18]; //   'x6' : 'dx7',
            yp[18] = a*(y_prev[16]-y[18]);//   'dx6'
        }
    }
};

// for DDE - variational equations for Jacobian with respect to the initial conditions have to be defined
class DdeJacFunc {
public:
    template <class C>
    //  void operator()(C yp[sysdim], C y[sysdim], int mode) {
    void operator()(vector<vector<C>> &Jp, vector<vector<C>> &J, vector<vector<C>> &J_prev, vector<C> &x,  vector<C> &x_prev) {
        // running example
        if (syschoice == 1)  // running example
        {
            Jp[0][0] = -J[0][0]*x_prev[0] - J_prev[0][0]*x[0];
        }
        else if (syschoice == 2) //  example 5.15
        {
            Jp[0][0] = J_prev[0][0] + J[1][0];
            Jp[0][1] = J_prev[0][1] + J[1][1];
            Jp[1][0] = J[0][0] - J_prev[0][0];
            Jp[1][1] = J[0][1] - J_prev[0][1];
        }
        else if (syschoice == 3) // Xue et al. 2017 ex 3
        {
            Jp[0][0] = 1.4*J[2][0] - 0.9*J_prev[0][0];
            Jp[0][1] = 1.4*J[2][1] - 0.9*J_prev[0][1];
            Jp[0][2] = 1.4*J[2][2] - 0.9*J_prev[0][2];
            Jp[0][3] = 1.4*J[2][3] - 0.9*J_prev[0][3];
            Jp[0][4] = 1.4*J[2][4] - 0.9*J_prev[0][4];
            Jp[0][5] = 1.4*J[2][5] - 0.9*J_prev[0][5];
            Jp[0][6] = 1.4*J[2][6] - 0.9*J_prev[0][6];
            
            Jp[1][0] = 2.5*J[4][0] - 1.5*J[1][0];
            Jp[1][1] = 2.5*J[4][1] - 1.5*J[1][1];
            Jp[1][2] = 2.5*J[4][2] - 1.5*J[1][2];
            Jp[1][3] = 2.5*J[4][3] - 1.5*J[1][3];
            Jp[1][4] = 2.5*J[4][4] - 1.5*J[1][4];
            Jp[1][5] = 2.5*J[4][5] - 1.5*J[1][5];
            Jp[1][6] = 2.5*J[4][6] - 1.5*J[1][6];
            
            Jp[2][0] = 0.6*J[6][0] - 0.8*J[2][0]*x[1] - 0.8*x[2]*J[1][0];
            Jp[2][1] = 0.6*J[6][1] - 0.8*J[2][1]*x[1] - 0.8*x[2]*J[1][1];
            Jp[2][2] = 0.6*J[6][2] - 0.8*J[2][2]*x[1] - 0.8*x[2]*J[1][2];
            Jp[2][3] = 0.6*J[6][3] - 0.8*J[2][3]*x[1] - 0.8*x[2]*J[1][3];
            Jp[2][4] = 0.6*J[6][4] - 0.8*J[2][4]*x[1] - 0.8*x[2]*J[1][4];
            Jp[2][5] = 0.6*J[6][5] - 0.8*J[2][5]*x[1] - 0.8*x[2]*J[1][5];
            Jp[2][6] = 0.6*J[6][6] - 0.8*J[2][6]*x[1] - 0.8*x[2]*J[1][6];
            
            Jp[3][0] = - 1.3*J[3][0]*x[2] - 1.3*x[3]*J[2][0];
            Jp[3][1] = - 1.3*J[3][1]*x[2] - 1.3*x[3]*J[2][1];
            Jp[3][2] = - 1.3*J[3][2]*x[2] - 1.3*x[3]*J[2][2];
            Jp[3][3] = - 1.3*J[3][3]*x[2] - 1.3*x[3]*J[2][3];
            Jp[3][4] = - 1.3*J[3][4]*x[2] - 1.3*x[3]*J[2][4];
            Jp[3][5] = - 1.3*J[3][5]*x[2] - 1.3*x[3]*J[2][5];
            Jp[3][6] = - 1.3*J[3][6]*x[2] - 1.3*x[3]*J[2][6];
            
            Jp[4][0] = 0.7*J[0][0] - J[3][0]*x[4] - x[3]*J[4][0];
            Jp[4][1] = 0.7*J[0][1] - J[3][1]*x[4] - x[3]*J[4][1];
            Jp[4][2] = 0.7*J[0][2] - J[3][2]*x[4] - x[3]*J[4][2];
            Jp[4][3] = 0.7*J[0][3] - J[3][3]*x[4] - x[3]*J[4][3];
            Jp[4][4] = 0.7*J[0][4] - J[3][4]*x[4] - x[3]*J[4][4];
            Jp[4][5] = 0.7*J[0][5] - J[3][5]*x[4] - x[3]*J[4][5];
            Jp[4][6] = 0.7*J[0][6] - J[3][6]*x[4] - x[3]*J[4][6];
            
            Jp[5][0] = 0.3*J[0][0] - 3.1*J[5][0];
            Jp[5][1] = 0.3*J[0][1] - 3.1*J[5][1];
            Jp[5][2] = 0.3*J[0][2] - 3.1*J[5][2];
            Jp[5][3] = 0.3*J[0][3] - 3.1*J[5][3];
            Jp[5][4] = 0.3*J[0][4] - 3.1*J[5][4];
            Jp[5][5] = 0.3*J[0][5] - 3.1*J[5][5];
            Jp[5][6] = 0.3*J[0][6] - 3.1*J[5][6];

            
            Jp[6][0] = 1.8*J[5][0] - 1.5*J[6][0]*x[1] - 1.5*x[6]*J[1][0];
            Jp[6][1] = 1.8*J[5][1] - 1.5*J[6][1]*x[1] - 1.5*x[6]*J[1][1];
            Jp[6][2] = 1.8*J[5][2] - 1.5*J[6][2]*x[1] - 1.5*x[6]*J[1][2];
            Jp[6][3] = 1.8*J[5][3] - 1.5*J[6][3]*x[1] - 1.5*x[6]*J[1][3];
            Jp[6][4] = 1.8*J[5][4] - 1.5*J[6][4]*x[1] - 1.5*x[6]*J[1][4];
            Jp[6][5] = 1.8*J[5][5] - 1.5*J[6][5]*x[1] - 1.5*x[6]*J[1][5];
            Jp[6][6] = 1.8*J[5][6] - 1.5*J[6][6]*x[1] - 1.5*x[6]*J[1][6];
        }
        else if (syschoice == 4) // Szczelina_1 2014
        {
            Jp[0][0] = -J[0][0] + J_prev[0][0]*x_prev[0];
        }
        else if (syschoice == 5) // Szczelina_2 2014
        {
            Jp[0][0] = -J[0][0] -3.2*J_prev[0][0] + 3.*J_prev[0][0]*x_prev[0]*x_prev[0];
        }
        else if (syschoice == 6)  // self-driving car
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[1][0] = - params[0]*J_prev[0][0] - params[1]*J_prev[1][0];
            Jp[1][1] = - params[0]*J_prev[0][1] - params[1]*J_prev[1][1];
        }
        else if (syschoice == 7)  // self-driving car with uncertain (but constant) coefficients
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[0][2] = J[1][2];
            Jp[0][3] = J[1][3];
            Jp[1][0] = -J[2][0]*x_prev[0] - x[2] *J_prev[0][0] + J[2][0] - J[3][0]*x_prev[1] - x[3]*J_prev[1][0];
            Jp[1][1] = -J[2][1]*x_prev[0] - x[2] *J_prev[0][1] + J[2][1]  - J[3][1]*x_prev[1] - x[3]*J_prev[1][1];
            Jp[1][2] = -J[2][2]*x_prev[0] - x[2] *J_prev[0][2] + J[2][2]  - J[3][2]*x_prev[1] - x[3]*J_prev[1][2];
            Jp[1][3] = -J[2][3]*x_prev[0] - x[2] *J_prev[0][3] + J[2][3]  - J[3][3]*x_prev[1] - x[3]*J_prev[1][3];
            for (int i=0 ; i<sysdim ; i++)
            {
                Jp[2][i] = 0;
                Jp[3][i] = 0;
            }
        }
        else if (syschoice == 8)  // self-driving car with uncertain (but constant) coefficients
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[0][2] = J[1][2];
            Jp[0][3] = J[1][3];
            Jp[1][0] = - inputs[2]*J_prev[0][0]   - inputs[3]*J_prev[1][0];
            Jp[1][1] = - inputs[2]*J_prev[0][1]   - inputs[3]*J_prev[1][1];
            Jp[1][2] = -x_prev[0] - inputs[2]*J_prev[0][2] + 1 - inputs[3]*J_prev[1][2];
            Jp[1][3] = - inputs[2]*J_prev[0][3] - x_prev[1] - inputs[3]*J_prev[1][3];
        }
        else if (syschoice == 9) // Ex 4 of Zou CAV 2015
        {
            Jp[0][0] = -3*x_prev[0]*x_prev[0]*J_prev[0][0];
        }
        else if (syschoice == 10) // platoon
        {
            double a = 2.5;
            //       yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
            Jp[0][0] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][0]/30;
            Jp[0][1] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][1]/30;
            Jp[0][2] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][2]/30;
            Jp[0][3] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][3]/30;
            Jp[0][4] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][4]/30;
            Jp[0][5] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][5]/30;
            Jp[0][6] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][6]/30;
            Jp[0][7] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][7]/30;
            Jp[0][8] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][8]/30;
            //      yp[1] = y[2]; //  'x2' : 'dx2',
            Jp[1][0] = J[2][0];
            Jp[1][1] = J[2][1];
            Jp[1][2] = J[2][2];
            Jp[1][3] = J[2][3];
            Jp[1][4] = J[2][4];
            Jp[1][5] = J[2][5];
            Jp[1][6] = J[2][6];
            Jp[1][7] = J[2][7];
            Jp[1][8] = J[2][8];
            //      yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6. - y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
            Jp[2][0] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][0]/30 - J_prev[2][0]);
            Jp[2][1] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][1]/30 - J_prev[2][1]);
            Jp[2][2] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][2]/30 - J_prev[2][2]);
            Jp[2][3] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][3]/30 - J_prev[2][3]);
            Jp[2][4] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][4]/30 - J_prev[2][4]);
            Jp[2][5] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][5]/30 - J_prev[2][5]);
            Jp[2][6] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][6]/30 - J_prev[2][6]);
            Jp[2][7] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][7]/30 - J_prev[2][7]);
            Jp[2][8] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][8]/30 - J_prev[2][8]);
             //      yp[3] = y[4]; //  'x3' : 'dx3',
            Jp[3][0] = J[4][0];
            Jp[3][1] = J[4][1];
            Jp[3][2] = J[4][2];
            Jp[3][3] = J[4][3];
            Jp[3][4] = J[4][4];
            Jp[3][5] = J[4][5];
            Jp[3][6] = J[4][6];
            Jp[3][7] = J[4][7];
            Jp[3][8] = J[4][8];
            //     yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
            Jp[4][0] = a*(J_prev[2][0]-J_prev[4][0]);
            Jp[4][1] = a*(J_prev[2][1]-J_prev[4][1]);
            Jp[4][2] = a*(J_prev[2][2]-J_prev[4][2]);
            Jp[4][3] = a*(J_prev[2][3]-J_prev[4][3]);
            Jp[4][4] = a*(J_prev[2][4]-J_prev[4][4]);
            Jp[4][5] = a*(J_prev[2][5]-J_prev[4][5]);
            Jp[4][6] = a*(J_prev[2][6]-J_prev[4][6]);
            Jp[4][7] = a*(J_prev[2][7]-J_prev[4][7]);
            Jp[4][8] = a*(J_prev[2][8]-J_prev[4][8]);
             //     yp[5] = y[6]; //   'x4' : 'dx4',
            Jp[5][0] = J[6][0];
            Jp[5][1] = J[6][1];
            Jp[5][2] = J[6][2];
            Jp[5][3] = J[6][3];
            Jp[5][4] = J[6][4];
            Jp[5][5] = J[6][5];
            Jp[5][6] = J[6][6];
            Jp[5][7] = J[6][7];
            Jp[5][8] = J[6][8];
            //     yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
            Jp[6][0] = a*(J_prev[4][0]-J_prev[6][0]);
            Jp[6][1] = a*(J_prev[4][1]-J_prev[6][1]);
            Jp[6][2] = a*(J_prev[4][2]-J_prev[6][2]);
            Jp[6][3] = a*(J_prev[4][3]-J_prev[6][3]);
            Jp[6][4] = a*(J_prev[4][4]-J_prev[6][4]);
            Jp[6][5] = a*(J_prev[4][5]-J_prev[6][5]);
            Jp[6][6] = a*(J_prev[4][6]-J_prev[6][6]);
            Jp[6][7] = a*(J_prev[4][7]-J_prev[6][7]);
            Jp[6][8] = a*(J_prev[4][8]-J_prev[6][8]);
             //     yp[7] = y[8];//  'x5' : 'dx5',
            Jp[7][0] = J[8][0];
            Jp[7][1] = J[8][1];
            Jp[7][2] = J[8][2];
            Jp[7][3] = J[8][3];
            Jp[7][4] = J[8][4];
            Jp[7][5] = J[8][5];
            Jp[7][6] = J[8][6];
            Jp[7][7] = J[8][7];
            Jp[7][8] = J[8][8];
             //      yp[8] = a*(y_prev[6]-y[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
            Jp[8][0] = a*(J_prev[6][0]-J[8][0]);
            Jp[8][1] = a*(J_prev[6][1]-J[8][1]);
            Jp[8][2] = a*(J_prev[6][2]-J[8][2]);
            Jp[8][3] = a*(J_prev[6][3]-J[8][3]);
            Jp[8][4] = a*(J_prev[6][4]-J[8][4]);
            Jp[8][5] = a*(J_prev[6][5]-J[8][5]);
            Jp[8][6] = a*(J_prev[6][6]-J[8][6]);
            Jp[8][7] = a*(J_prev[6][7]-J[8][7]);
            Jp[8][8] = a*(J_prev[6][8]-J[8][8]);
        }
        else if (syschoice == 11) // platoon of 10 vehicles
        {
            double a = 2.5;
           
            for (int j=0; j<jacdim ; j++)
            {
                 //       yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
                Jp[0][j] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][j]/30;
             //      yp[1] = y[2]; //  'x2' : 'dx2',
                Jp[1][j] = J[2][j];
            //      yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6. - y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
                Jp[2][j] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][j]/30 - J_prev[2][j]);
            //      yp[3] = y[4]; //  'x3' : 'dx3',
                Jp[3][j] = J[4][j];
            //     yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
                Jp[4][j] = a*(J_prev[2][j]-J_prev[4][j]);
            //     yp[5] = y[6]; //   'x4' : 'dx4',
                Jp[5][j] = J[6][j];
            //     yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
                Jp[6][j] = a*(J_prev[4][j]-J_prev[6][j]);
            //     yp[7] = y[8];//  'x5' : 'dx5',
                Jp[7][j] = J[8][j];
            //      yp[8] = a*(y_prev[6]-y[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
                Jp[8][j] = a*(J_prev[6][j]-J_prev[8][j]);
                //
                Jp[9][j] = J[10][j];
                Jp[10][j] = a*(J_prev[8][j]-J_prev[10][j]);
                Jp[11][j] = J[12][j];
                Jp[12][j] = a*(J_prev[10][j]-J_prev[12][j]);
                Jp[13][j] = J[14][j];
                Jp[14][j] = a*(J_prev[12][j]-J_prev[14][j]);
                Jp[15][j] = J[16][j];
                Jp[16][j] = a*(J_prev[14][j]-J_prev[16][j]);
                Jp[17][j] = J[18][j];
                Jp[18][j] = a*(J_prev[16][j]-J[18][j]);
            }
        }
    }
};






#endif
