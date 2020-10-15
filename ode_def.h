/* ============================================================================
 File   : ode_def.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 The place to declare the systems of ODEs or DDEs on which he/she wants to perform reachability
 Setup of dimension of the system (sysdim), initial conditions and parameters is in ode_def.cpp
 ============================================================================ */
#ifndef ODE_DEF_H
#define ODE_DEF_H

#include "utils.h"
#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"

using namespace std;

#define innerapprox 1

extern int systype;    // systype = 0 (ODE) or 1 (DDE) -- initialized in main.cpp / command line
extern int syschoice;  // to choose among the predefined systems of ODE or DDE -- initialized in main.cpp / command line


extern int sysdim; // dimension of the system of ODE/DDE to analyze
extern int inputsdim; // dimension of the uncertain inputs and parameters of the system
extern int fullinputsdim; // full dimension of the uncertain inputs and parameters of the system: taking into account variable inputs
extern int jacdim; // Jacobian will be dimension jacdim = sysdim + fullinputsdim
extern int sysdim_params;  // dimension of the vector of parameters params that do not appear in Jacobian

extern double t_end; // ending time of integration
extern double t_begin; // starting time of initialization

extern vector<AAF> params;      // params of the ODE (nondeterministic disturbances)

extern vector<AAF> initial_values; // uncertain initial conditions
extern vector<AAF> center_initial_values;

extern vector<AAF> inputs;   // uncertain inputs and parameters
extern vector<AAF> fullinputs; // uncertain inputs and parameters
extern vector<int> nb_inputs; // piecewise constant input changes value every t_end/nb_inputs[i] seconds

extern vector<AAF> center_fullinputs;
extern vector<int> index_param;
extern vector<int> index_param_inv;
extern vector<interval> eps;

extern vector<vector<interval>> Jac_param_inputs; // for inputs defined as g(x1,...xn): we give the jacobian

// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern int component_to_subdiv, component_to_subdiv2;

extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_joint_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
extern vector<double> t_print; // times where results are stored
extern int current_subdiv;
extern int current_iteration;

// for robust inner-approximations
extern int uncontrolled;  // number of uncontrolled parameters (forall params)
extern int controlled;  // number of controlled parameters (forall params)
extern vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
//extern vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust inner-approx)
//extern int variable;  // number of non constant parameters
//extern vector<bool> is_variable; // for each parameter, constant or variable

extern bool refined_mean_value;

extern bool print_debug;


void define_system_dim(int argc, char* argv[]);  // define the dimensions of your system (ODE or DDE)



// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &params_inputs, vector<AAF> &param_inputs_center, vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for DDEs : functions that initialize the DDE on first time period
vector<T<AAF>> Initfunc(const T<AAF>& t, vector<AAF> &beta_initial, vector<AAF> &beta_inputs);
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta_initial, vector<T<F<AAF>>> &beta_inputs);

// defining analytical solution if any for comparison
void AnalyticalSol(int current_iteration, double d0);

// reading sysdim, jacdim, etc
void readfromfile_system_dim(const char * params_filename, int &sysdim, int &jacdim, int &sysdim_params, int &nb_subdiv_init);

// d0 and t_begin and nb_subdiv are for DDEs only, rest are common to ODE and DDE
void read_parameters(const char * params_filename, double &tau, double &t_end, double &d0, double &t_begin, int &order, int &nb_subdiv);

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(int argc, char* argv[], double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

void init_utils_inputs(double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> initial_values_save, vector<AAF> inputs_save, int param_to_subdivide);


// define here  your ODE system yp = \dot y = f(y)
class OdeFunc {
public:
    template <class C>
      void operator()(vector<C> &yp, vector<C> param_inputs, vector<C> y) {
          
          
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
              yp[0] = - g*(y[1])-rho*y[0]*y[0]*a*d/(2.0*param_inputs[0]); // velocity v
           yp[1] = - g*(1-y[1]*y[1]/2)/y[0]; // angle gamma with respect to the x axis
           yp[2] = y[0]*(1-y[1]*y[1]/2); // position x
           yp[3] = y[0]*(y[1]); // position y
            //  yp[4] = 0;
          }
          /******************* end ballistic ****************************/
          else if (syschoice == 5)  // self-driving car
          {
              yp[0] = y[1];
              yp[1] = -params[0] *(y[0] - 1.0) - params[1]*y[1];   // pr = 1 is the reference position
          }
          else if (syschoice == 6)  // self-driving car with constant parameters
          {
              yp[0] = y[1];
              yp[1] = -param_inputs[0] *(y[0] - 1.0) - param_inputs[1]*y[1];   // pr = 1 is the reference position
          }
          else if (syschoice == 7)  // self-driving car with time-varying parameters
          {
              yp[0] = y[1];
              yp[1] = -y[2] *(y[0] - 1.0) - y[3]*y[1];   // pr = 1 is the reference position
              yp[2] = interval(-2.,2.); //AAF(interval(-2.,2.)) + y[4];  // parameter Kp
              yp[3] = interval(-2.,2.); // 0;  // parameter Kd
          }
          else if (syschoice == 8)
          {
              yp[0] = - y[0]*y[0]*y[0];
          }
          else if (syschoice == 9)  // acrobatic quadcopter
          {
              double Cv = 0.25;     // transitional drag
              double Cphi = 0.02255; // rotational drag
              double g = 9.81;
              double m = 1.25;      // mass
              double l = 0.5;       // length of quadrotor"s center to an edge
              double Iyy = 0.03;    // moment of inertia
              yp[0] = -y[1];     // px
              yp[1] = -(-Cv/m*y[1] - (param_inputs[0]+param_inputs[1])/m*sin(y[4]));     // vx
              yp[2] = -y[3];     // py
              yp[3] = -(-(g+Cv*y[2]/m) + (param_inputs[0]+param_inputs[1])/m*cos(y[4]));     // vy
              yp[4] = -y[5];     // phi
              yp[5] = -(-Cphi/Iyy*y[5] + l/Iyy*(param_inputs[1]-param_inputs[0]));     // omega
          }
          else if (syschoice == 99)  // acrobatic quadcopter with m et Iyy as disturbances
          {
              double Cv = 0.25;     // transitional drag
              double Cphi = 0.02255; // rotational drag
              double g = 9.81;
              // double m = 1.25;      // mass : param_inputs[2]
              double l = 0.5;       // length of quadrotor"s center to an edge
              // double Iyy = 0.03;    // moment of inertia : param_inputs[3]
              yp[0] = -y[1];     // px
              yp[1] = -(-Cv/param_inputs[2]*y[1] - (param_inputs[0]+param_inputs[1])/param_inputs[2]*sin(y[4]));     // vx
              yp[2] = -y[3];     // py
              yp[3] = -(-(g+Cv*y[2]/param_inputs[2]) + (param_inputs[0]+param_inputs[1])/param_inputs[2]*cos(y[4]));     // vy
              yp[4] = -y[5];     // phi
              yp[5] = -(-Cphi/param_inputs[3]*y[5] + l/param_inputs[3]*(param_inputs[1]-param_inputs[0]));     // omega
          }
          else if (syschoice == 10)  // 10-D near-hover quadrotor
          {
              double n0 = 10;
              double d0 = 10;
              double d1 = 8;
              double kT = 0.91;
              double m = 1.3;
              double g = 9.81;
              // param_inputs[0] .. param_inputs[2] : disturbances dx, dy, dz
              // param_inputs[3] .. param_inputs[5] : control inputs Sx, Sy, Tz
              yp[0] = -(y[1] + param_inputs[0]);       // px' = vx + dx
              yp[1] = -(g*tan(y[2]));                  // vx' = g.tan(thetax)
              yp[2] = -(-d1*y[2] +y[3]);               // thetax' = -d1.thetax + omegax
              yp[3] = -(-d0*y[2] + n0*param_inputs[3]); // omegax' = -d0.thetax + no.Sx
              yp[4] = -(y[5] + param_inputs[1]);        // py' = vy + dy
              yp[5] = -(g*tan(y[6]));                   // vy' = g.tan(thetay)
              yp[6] = -(-d1*y[6] +y[7]);                // thetay' = -d1.thetay + omegay
              yp[7] = -(-d0*y[6] + n0*param_inputs[4]);  // omegay' = -d0.thetay + no.Sy
              yp[8] = -(y[9] + param_inputs[2]);        // pz' = vz + dz
              yp[9] = -(kT*param_inputs[5] - g);        // vz' = kT.Tz - g
          }
          else if (syschoice == 11)  // Dubbins vehicle
          {
              double v = 5;         // constant velocity
              // param_inputs[0] .. param_inputs[2] : disturbances b1, b2, b3
              // param_inputs[3]  : control input a : angular control
              yp[0] = -(v*cos(y[2]) + param_inputs[0]);        // px' = v.cos(theta) + b1
              yp[1] = -(v*sin(y[2]) + param_inputs[1]);        // py' = v.sin(theta) + b2
              yp[2] = -(param_inputs[3] + param_inputs[2]);    // theta' = a + b3
          }
          else if (syschoice == 12)  // academic example to investigate time-varying parameters
          {
              yp[0] = (param_inputs[0] + param_inputs[1]*y[1])*y[0];
              yp[1] = 1;
          }
          else if (syschoice == 13)  // Laub-Loomis Benchmark [Arch 2019]
          {
              yp[0] = 1.4*y[2] - 0.9*y[0];
              yp[1] = 2.5*y[4] - 1.5*y[1];
              yp[2] = 0.6*y[6] - 0.8*y[1]*y[2];
              yp[3] = 2 - 1.3*y[2]*y[3];
              yp[4] = 0.7*y[0] - y[3]*y[4];
              yp[5] = 0.3*y[0] - 3.1*y[5];
              yp[6] = 1.8*y[5]-1.5*y[1]*y[6];
          }
          else if (syschoice == 14) // Van der Pol oscillator [Arch 2019]
          {
              double mu = 1.;
              yp[0] = y[1];
              yp[1] = mu*(1-y[0]*y[0])*y[1]-y[0];
          }
          else if (syschoice == 15) // Van der Pol oscillator [Arch 2018 and Sparse Polynomial zonotopes]
          {
              yp[0] = y[1];
              yp[1] = (1-y[0]*y[0])*y[1]-y[0];
          }
          else if(syschoice == 17) // quadrotor model [Arch 2019]
          {
              double g = 9.81;
              double R = 0.1;
              double l = 0.5;
              double Mmotor = 0.5;
              double M = 1;
              double m = M + 4*Mmotor;
              double Jx = 2./5*M*R*R + 2*l*l*Mmotor;
              double Jy = Jx;
              double Jz = 2./5*M*R*R + 4*l*l*Mmotor;
              double u1 = 1.;
              double u2 = 0;
              double u3 = 0;
              auto F = m*g - 10*(y[2]-u1)+3*y[5]; // height control
              auto tau_phi = -(y[6]-u2)-y[9];  // roll control
              auto tau_theta = -(y[7]-u3)-y[10]; // pitch control
              auto tau_psi = 0; // heading uncontrolled
              yp[0] = cos(y[7])*cos(y[8])*y[3] + (sin(y[6])*sin(y[7])*cos(y[8])-cos(y[6]*sin(y[8])))*y[4] + (cos(y[6])*sin(y[7])*cos(y[8])+sin(y[6])*sin(y[8]))*y[5];
              yp[1] = cos(y[7])*sin(y[8])*y[3] + (sin(y[6])*sin(y[7])*cos(y[8])+cos(y[6]*sin(y[8])))*y[4] + (cos(y[6])*sin(y[7])*sin(y[8])-sin(y[6])*cos(y[8]))*y[5];
              yp[2] = sin(y[7])*y[3] - sin(y[6])*cos(y[7])*y[4] - cos(y[6])*cos(y[7])*y[5];
              yp[3] = y[11]*y[4] - y[10]*y[5] - g*sin(y[7]);
              yp[4] = y[9]*y[5] - y[11]*y[3] + g*cos(y[7])*sin(y[6]);
              yp[5] = y[10]*y[3] - y[9]*y[4] + g*cos(y[7])*cos(y[6]) - F/m;
              yp[6] = y[9] + sin(y[6])*tan(y[7])*y[10] + cos(y[6])*tan(y[7])*y[11];
              yp[7] = cos(y[6])*y[10] - sin(y[6])*y[11];
              yp[8] = sin(y[6])/cos(y[7])*y[10] + cos(y[6])/cos(y[7])*y[11];
              yp[9] = (Jy-Jz)/Jx*y[10]*y[11] + tau_phi/Jx;
              yp[10] = (Jz-Jx)/Jy*y[9]*y[11] + tau_theta/Jy;
              yp[11] = (Jx-Jy)/Jz*y[9]*y[10] + tau_psi/Jz;
          }
         else if(syschoice == 18) // HSCC 2019 paper crazyflie example
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
         else if (syschoice == 181) // HSCC 2019 paper crazyflie example now controlled by neural network
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
             
             auto sqp0 = param_inputs[0]*param_inputs[0]; // exp(2.0*param_inputs[0]);
	     auto sqp1 = param_inputs[1]*param_inputs[1]; // exp(2.0*param_inputs[1]);
	     auto sqp2 = param_inputs[2]*param_inputs[2]; // exp(2.0*param_inputs[2]);
             
             auto cmd_r = 400*param_inputs[0]*(10+sqp0)*(60+sqp0)/(600+sqp0*(270+sqp0*(11+sqp0/24))); // 800.*(expp0-1.0)/(expp0+1.0); // cmd_phi = 800*tanh(param_inputs[0]) // y[6]*Ki_rr + err_p*Kp_rr;
	     auto cmd_p = 400*param_inputs[1]*(10+sqp1)*(60+sqp1)/(600+sqp1*(270+sqp1*(11+sqp1/24))); // 800.0*(expp1-1.0)/(expp1+1.0); // cmd_theta // y[7]*Ki_pr + err_q*Kp_pr;
	     auto cmd_y = 1000*param_inputs[2]*(10+sqp2)*(60+sqp2)/(600+sqp2*(270+sqp2*(11+sqp2/24))); // 3000.0*(expp2-1.0)/(expp2+1.0); // cmd_psi // y[8]*Ki_yr + err_r*Kp_yr;
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
         else if (syschoice == 182) // HSCC 2019 paper crazyflie example but with a different altitude controller - and now controlled by neural network
         {
             static const double p_sp = 0.0;  // angular speed of 1 degree / sec
             static const double q_sp = 0.0;
             static const double r_sp = 0.2; // prefer yaw command
             
             static const double z_sp = 0.0; // just stabilize at 0
             
             static const double C1 = 0.04076521;
             static const double C2 = 380.8359;
             static const double d = 0.046/sqrt(2.0);
             
             static const double Ixx = 1.657171e-5;
             static const double Iyy = 1.6655602e-5;
             static const double Izz = 2.9261652e-5;
             
             static const AAF Ct = 1.285e-8;//interval(1.28e-8 , 1.29e-8);
             static const AAF Cd = 7.645e-11;//interval(7.64e-11 , 7.65e-11);
             
             //static const double Kp_Z = 2.0;
             //static const double Kp_VZ = 25.0;
             //static const double Ki_VZ = 15.0;
             static const double Kp_Z = 3000.0;
             static const double Kp_VZ = 500.0;
             static const double Ki_VZ = 300.0;
             
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
             //auto thrust_Raw     = Kp_VZ * (yp[13]) + Ki_VZ * y[13];
             
             auto thrust = Kp_VZ * (yp[13]) + Ki_VZ * y[13] + Kp_VZ * y[11] + 48500.0; //1000.0*thrust_Raw + 36000;
             
             auto err_p = p_sp - y[3];
             auto err_q = q_sp - y[4];
             auto err_r = r_sp - y[5];
             
             auto sqp0 = param_inputs[0]*param_inputs[0]; // exp(2.0*param_inputs[0]);
	     auto sqp1 = param_inputs[1]*param_inputs[1]; // exp(2.0*param_inputs[1]);
	     auto sqp2 = param_inputs[2]*param_inputs[2]; // exp(2.0*param_inputs[2]);
             
             auto cmd_r = 400*param_inputs[0]*(10+sqp0)*(60+sqp0)/(600+sqp0*(270+sqp0*(11+sqp0/24))); // 800.*(expp0-1.0)/(expp0+1.0); // cmd_phi = 800*tanh(param_inputs[0]) // y[6]*Ki_rr + err_p*Kp_rr;
	     auto cmd_p = 400*param_inputs[1]*(10+sqp1)*(60+sqp1)/(600+sqp1*(270+sqp1*(11+sqp1/24))); // 800.0*(expp1-1.0)/(expp1+1.0); // cmd_theta // y[7]*Ki_pr + err_q*Kp_pr;
	     auto cmd_y = 1000*param_inputs[2]*(10+sqp2)*(60+sqp2)/(600+sqp2*(270+sqp2*(11+sqp2/24))); // 3000.0*(expp2-1.0)/(expp2+1.0); // cmd_psi // y[8]*Ki_yr + err_r*Kp_yr;
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
         else if (syschoice == 19) {  // academic example, time-varying (piecewise constant) parameters
             yp[0] = 2 + 2*param_inputs[0]*(1-y[1]) + y[1];
             yp[1] = 1; // time
         }
         
         else if (syschoice == 21) {  // academic example, time-varying (piecewise constant) parameters
             yp[0] = 2 + 2*param_inputs[0]*(1-y[1]) + y[1];
             yp[1] = 1; // time
         }
         else if (syschoice == 22) {  // academic example, time-varying (piecewise constant) parameters
             yp[0] = 2 + 2*(param_inputs[0]+param_inputs[0]*param_inputs[0])*(1-y[1]) + y[1];
             yp[1] = 1; // time
         }
         else if (syschoice == 23) {   // pursuer-evader example Mitchell
             double v = 5.0;  // linear velocity
             // param_inputs[0] : angular velocity a of the evader (control)
             // param_inputs[1] : angular velocity b of the pursuer (disturbance)
             yp[0] = -v + v*cos(y[2]) + param_inputs[0]*y[1];
             yp[1] = v*sin(y[2]) - param_inputs[0]*y[1];
             yp[2] = param_inputs[1] - param_inputs[0];
         }
         else if (syschoice == 24) { // [Franzle et al.]
             yp[0] = -(-0.5*y[0] - (0.5+param_inputs[0])*y[1] + 0.5);
             yp[1] = - (-0.5*y[1] + 1.0);
         }
         else if (syschoice == 25) { // [Franzle et al.] reversed time van der pol oscillator with uncertainty
             yp[0] = y[1];
             yp[1] = - (0.4*y[0] + 5.0*(y[0]*y[0] - (param_inputs[0] + 0.2))*y[1]);
         }
         else if (syschoice == 26) { // [Franzle et al.] 7-d biological system
             yp[0] = -(-0.4*y[0] + 5.0*y[2]*y[3] + param_inputs[0]);
             yp[1] = -(0.4*y[0] - y[1]);
             yp[2] = -(y[1] - 5.0*y[2]*y[3]);
             yp[3] = -(5*y[4]*y[5] - 5.0*y[2]*y[3]);
             yp[4] = -(-5*y[4]*y[5] + 5.0*y[2]*y[3]);
             yp[5] = -(0.5*y[6] - 5.0*y[4]*y[5]);
             yp[6] = -(-0.5*y[6] + 5.0*y[4]*y[5]);
         }
         else if (syschoice == 27) { // [Franzle et al.] 7-d biological system but with sharing
             C y2y3 = y[2]*y[3];
             C y4y5 = y[4]*y[5];
             yp[0] = -(-0.4*y[0] + 5.0*y2y3 + param_inputs[0]);
             yp[1] = -(0.4*y[0] - y[1]);
             yp[2] = -(y[1] - 5.0*y2y3);
             yp[3] = -(5*y4y5 - 5.0*y2y3);
             yp[4] = -(-5*y4y5 + 5.0*y2y3);
             yp[5] = -(0.5*y[6] - 5.0*y4y5);
             yp[6] = -(-0.5*y[6] + 5.0*y4y5);
         }
         else if (syschoice == 28) { // [Franzle et al.] 7-d biological system without disturbance
             C y2y3 = y[2]*y[3];
             C y4y5 = y[4]*y[5];
             yp[0] = -(-0.4*y[0] + 5.0*y2y3);
             yp[1] = -(0.4*y[0] - y[1]);
             yp[2] = -(y[1] - 5.0*y2y3);
             yp[3] = -(5*y4y5 - 5.0*y2y3);
             yp[4] = -(-5*y4y5 + 5.0*y2y3);
             yp[5] = -(0.5*y[6] - 5.0*y4y5);
             yp[6] = -(-0.5*y[6] + 5.0*y4y5);
         }
         else if (syschoice == 29) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
             yp[0] = y[3]*cos(y[2]);
             yp[1] = y[3]*sin(y[2]);
             yp[2] = param_inputs[0];
             yp[3] = param_inputs[1] + interval(-0.0001,0.0001);
         }
        else if (syschoice == 30) { // EX_1 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            yp[0] = y[1]-y[0]*y[0]*y[0];
            yp[1] = param_inputs[0];
        }
        else if (syschoice == 31) { // Quadcopter Mikhail Bessa avec 3 composantes de plus
            // y[0]=z ; y[1]=u ; y[2]=v ; y[3]=w ; y[4]=phi ; y[5]=theta ;
            // y[6]=psi ; y[7]=p ; y[8]=q ; y[9]=r ; y[10]=zI ;
            // param_inputs[0]=cmd_phi ; param_inputs[1]=cmd_theta ;
            // param_inputs[2]=cmd_psi ;
            yp[0] = (-sin(y[5]))*y[1] + cos(y[5])*sin(y[4])*y[2] + cos(y[5])*cos(y[4])*y[3]; // yp[0] = z
            yp[1] = y[9]*y[2] - y[8]*y[3] + sin(y[5])*9.81; // yp[1] = u
            yp[2] = -y[9]*y[1] + y[7]*y[3] - cos(y[5])*sin(y[4])*9.81; // yp[2] = v
            yp[3] = y[8]*y[1] - y[7]*y[2] - 9.81*cos(y[5])*cos(y[4])+(7.62648576804e-10)*param_inputs[1]*param_inputs[1] + (7.62648576804e-10)*param_inputs[0]*param_inputs[0] + (3.05059430721e-9)*param_inputs[2]*param_inputs[2] + 7.62648576804*y[0]*y[0] + 7.62648576804*y[3]*y[0] - 4.57589146082*y[0]*y[10] + 1.90662144201*y[3]*y[3] - 2.28794573041*y[3]*y[10] + 0.686383719125*y[10]*y[10] - 29.0850309334*y[0] - 14.5425154667*y[3] + 8.72550928*y[10] + 27.7303023346; // yp[3] = w
            yp[4] = y[7] + sin(y[5])/cos(y[5])*(cos(y[4])*y[9] + sin(y[4])*y[8]); // yp[4] = phi
            yp[5] = cos(y[4])*y[8] - sin(y[4])*y[9]; // yp[5] = theta
            yp[6] = cos(y[4]) / cos(y[5]) * y[9] + sin(y[4]) / cos(y[5]) * y[8]; // yp[6] = psi
            yp[7] = -0.760697* y[8] * y[9] + (-8.38277868313e-3)*y[0]*param_inputs[0]  - (4.19138934156e-3)*y[3]*param_inputs[0] + (2.51483360494e-3)*y[10]*param_inputs[0] + 0.0159846477606*param_inputs[0] - (1.67655573663e-7)*param_inputs[1]*param_inputs[2]; // yp[7] = p
            yp[8] = 0.761902* y[7] * y[9] + (-8.34055576802e-3)*y[0]*param_inputs[1] - (4.17027788401e-3)*y[3]*param_inputs[1] + (2.50216673041e-3)*y[10]*param_inputs[1] + 0.0159041352657*param_inputs[1] - (1.66811115361e-7)*param_inputs[0]*param_inputs[2]; // yp[8] = q
            yp[9] = -0.002867* y[7] * y[8] + (-1.73667282186e-3)*y[0]*param_inputs[2] - (8.68336410928e-4)*y[3]*param_inputs[2] + (5.21001846557e-4)*y[10]*param_inputs[2] + (3.31156343047e-3)*param_inputs[2] - (8.68336410931e-9)*param_inputs[1]*param_inputs[0]; // yp[9] = r
            yp[10] = 2*(1-y[0])-y[3]; // yp[10] = zI
            yp[11] = 0;
            yp[12] = 0;
            yp[13] = 0;
        }
        else if (syschoice == 311) { // Quadcopter Mikhail Bessa
            // y[0]=z ; y[1]=u ; y[2]=v ; y[3]=w ; y[4]=phi ; y[5]=theta ;
            // y[6]=psi ; y[7]=p ; y[8]=q ; y[9]=r ; y[10]=zI ;
            // param_inputs[0]=cmd_phi ; param_inputs[1]=cmd_theta ;
            // param_inputs[2]=cmd_psi ;
            yp[0] = (-sin(y[5]))*y[1] + cos(y[5])*sin(y[4])*y[2] + cos(y[5])*cos(y[4])*y[3]; // yp[0] = z
            yp[1] = y[9]*y[2] - y[8]*y[3] + sin(y[5])*9.81; // yp[1] = u
            yp[2] = -y[9]*y[1] + y[7]*y[3] - cos(y[5])*sin(y[4])*9.81; // yp[2] = v
            yp[3] = y[8]*y[1] - y[7]*y[2] - 9.81*cos(y[5])*cos(y[4])+(7.62648576804e-10)*param_inputs[1]*param_inputs[1] + (7.62648576804e-10)*param_inputs[0]*param_inputs[0] + (3.05059430721e-9)*param_inputs[2]*param_inputs[2] + 7.62648576804*y[0]*y[0] + 7.62648576804*y[3]*y[0] - 4.57589146082*y[0]*y[10] + 1.90662144201*y[3]*y[3] - 2.28794573041*y[3]*y[10] + 0.686383719125*y[10]*y[10] - 29.0850309334*y[0] - 14.5425154667*y[3] + 8.72550928*y[10] + 27.7303023346; // yp[3] = w
            yp[4] = y[7] + sin(y[5])/cos(y[5])*(cos(y[4])*y[9] + sin(y[4])*y[8]); // yp[4] = phi
            yp[5] = cos(y[4])*y[8] - sin(y[4])*y[9]; // yp[5] = theta
            yp[6] = cos(y[4]) / cos(y[5]) * y[9] + sin(y[4]) / cos(y[5]) * y[8]; // yp[6] = psi
            yp[7] = -0.760697* y[8] * y[9] + (-8.38277868313e-3)*y[0]*param_inputs[0]  - (4.19138934156e-3)*y[3]*param_inputs[0] + (2.51483360494e-3)*y[10]*param_inputs[0] + 0.0159846477606*param_inputs[0] - (1.67655573663e-7)*param_inputs[1]*param_inputs[2]; // yp[7] = p
            yp[8] = 0.761902* y[7] * y[9] + (-8.34055576802e-3)*y[0]*param_inputs[1] - (4.17027788401e-3)*y[3]*param_inputs[1] + (2.50216673041e-3)*y[10]*param_inputs[1] + 0.0159041352657*param_inputs[1] - (1.66811115361e-7)*param_inputs[0]*param_inputs[2]; // yp[8] = q
            yp[9] = -0.002867* y[7] * y[8] + (-1.73667282186e-3)*y[0]*param_inputs[2] - (8.68336410928e-4)*y[3]*param_inputs[2] + (5.21001846557e-4)*y[10]*param_inputs[2] + (3.31156343047e-3)*param_inputs[2] - (8.68336410931e-9)*param_inputs[1]*param_inputs[0]; // yp[9] = r
            yp[10] = 2*(1-y[0])-y[3]; // yp[10] = zI
        }
          else if (syschoice == 32) { // EX_2 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[1];
              yp[1] = param_inputs[0]*y[1]*y[1] - y[0];
          }
          else if (syschoice == 33) { // EX_3 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = -y[0]*(0.1 + (y[0]+y[1])*(y[0]+y[1]));
              yp[1] = (param_inputs[0]+y[0])*(0.1 + (y[0] + y[1])*(y[0]+y[1]));
          }
          else if (syschoice == 34) { // EX_4 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[1] + 0.5*y[2]*y[2];
              yp[1] = y[2];
              yp[2] = param_inputs[0];
          }
          else if (syschoice == 35) { // EX_5 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = -y[0] + y[1] - y[2];
              yp[1] = -y[0]*(y[2]+1) - y[1];
              yp[2] = -y[0] + param_inputs[0];
          }
          else if (syschoice == 36) { // EX_6 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = -y[0]*y[0]*y[0] + y[1];
              yp[1] = y[1]*y[1]*y[1] + y[2];
              yp[2] = param_inputs[0];
          }
          else if (syschoice == 37) { // EX_7 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[2]*y[2]*y[2] - y[1];
              yp[1] = y[2];
              yp[2] = param_inputs[0];
          }
          else if (syschoice == 38) { // EX_8 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[1];
              yp[1] = -9.8*y[2] + 1.6*y[2]*y[2]*y[2] + y[0]*y[3]*y[3];
              yp[2] = y[3];
              yp[3] = param_inputs[0];
          }
          else if (syschoice == 39) { // EX_9 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[1];
              yp[1] = -y[0] + 0.1*sin(y[2]);
              yp[2] = y[3];
              yp[3] = param_inputs[0];
          }
          else if (syschoice == 40) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
              yp[0] = y[3] * cos(y[2]);
              yp[1] = y[3] * sin(y[2]);
              yp[2] = param_inputs[1];
              yp[3] = param_inputs[0];
          }
          else if (syschoice == 41) { // essai sys couple
              yp[0] = y[1]+y[0];
              yp[1] = param_inputs[0];
          }
          else if(syschoice == 42) // HSCC 2019 paper crazyflie example+effets aerodynamiques
          {
              static const double p_sp = 1.0*M_PI/180.0;  // angular speed of 1 degree / sec
              static const double q_sp = 0.0;
              static const double r_sp = 0.0;
              
              static const double z_sp = 1.0;
              
              static const double C1 = 0.04076521;
              static const double C2 = 380.8359;
              static const double d = 0.046/sqrt(2.0);
              static const double h = 0.005;
              
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
              
              /* Cross-model for aerodynamical systems from Forster */
              /* 
              static const double K11 = -10.2506e-7; // e-8 and below?
              static const double K12 = -0.3177e-7;
              static const double K13 = -0.4332e-7;
              
              static const double K21 = -0.3177e-7;
              static const double K22 = -10.2506e-7;
              static const double K23 = -0.4332e-7;
              
              static const double K31 = -7.7050e-7;
              static const double K32 = -7.7050e-7;
              static const double K33 = -7.5530e-7;
	      */
	      
	      /* Simplified model */ 
	      static const double K11 = -9.1785e-8;
	      static const double K12 = 0;
	      static const double K13 = 0;
	      static const double K21 = 0;
	      static const double K22 = -9.1785e-8;
	      static const double K23 = 0;
	      static const double K31 = 0;
	      static const double K32 = 0;
	      static const double K33 = -10.311e-7;
	      
              /* Crazyflie trajectory tracking article*/
              static const AAF Ip_qr = (Iyy-Izz)/Ixx;//interval(-1.04880447793, -1.03580464787);
              static const AAF Iq_pr = (Izz-Ixx)/Iyy;//interval(1.03470095927, 1.04749270535);
              static const AAF Ir_pq = (Ixx-Iyy)/Izz;//interval(-0.0162919189567, -0.0120891632629);
              
              static const AAF Im_xx = 1.0/Ixx;//interval(71484.0524534, 71885.7226787);
              static const AAF Im_yy = 1.0/Iyy;//interval(69441.6509547, 69834.7034512);
              static const AAF Im_zz = 1.0/Izz;//interval(34492.4780616, 34712.0265858);
              
              static const double g = 9.8;
              static const double m = 0.028;
              
              auto cosRoll = cos(y[0]); // cos phi
              auto sinRoll = sin(y[0]); // sin phi
              
              auto cosPitch = cos(y[1]); // cos theta
              auto sinPitch = sin(y[1]); // sin theta
              
              auto cosYaw = cos(y[2]); // cos psi
              auto sinYaw = sin(y[2]); // sin psi
              
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
              
              auto cmd_r = y[6]*Ki_rr + err_p*Kp_rr; // cmd_phi
              auto cmd_p = y[7]*Ki_pr + err_q*Kp_pr; // cmd_theta
              auto cmd_y = y[8]*Ki_yr + err_r*Kp_yr; // cmd_psi
              
              /* computation of the PWMs */
              auto PWM1 = thrust-cmd_r/2-cmd_p/2-cmd_y;
              auto PWM2 = thrust-cmd_r/2+cmd_p/2+cmd_y;
              auto PWM3 = thrust+cmd_r/2+cmd_p/2-cmd_y;
              auto PWM4 = thrust+cmd_r/2-cmd_p/2+cmd_y;
              
              /* computaton of the angular velocities of the four rotors */
              auto Om1 = C1*PWM1+C2;
              auto Om2 = C1*PWM2+C2;
              auto Om3 = C1*PWM3+C2;
              auto Om4 = C1*PWM4+C2;
              
              /* R transpose, from inertial to body frame phi=Roll, theta=Pitch, psi=Yaw */
              /* first line */
              auto RT11 = cosYaw*cosPitch;
              auto RT12 = cosPitch*sinYaw;
              auto RT13 = -sinPitch;
              /* second line */
              auto RT21 = cosYaw*sinPitch*sinRoll-cosRoll*sinYaw;
              auto RT22 = cosYaw*cosRoll+sinYaw*sinPitch*sinRoll;
              auto RT23 = cosPitch*sinRoll;
              /* third line */
              auto RT31 = sinYaw*sinRoll+cosYaw*cosRoll*sinPitch;
              auto RT32 = cosRoll*sinYaw*sinPitch-cosYaw*sinPitch;
              auto RT33 = cosPitch*cosRoll;
              
              /* linear velocities of the four rotors */
	      auto u1 = y[9]+y[4]*h+y[5]*d*sqrt(2)/2;
	      auto u2 = y[9]+y[4]*h+y[5]*d*sqrt(2)/2;
	      auto u3 = y[9]+y[4]*h-y[5]*d*sqrt(2)/2;
	      auto u4 = y[9]+y[4]*h-y[5]*d*sqrt(2)/2;
  
	      auto v1 = y[10]-y[3]*h-y[5]*d*sqrt(2)/2;
	      auto v2 = y[10]-y[3]*h-y[5]*d*sqrt(2)/2;
	      auto v3 = y[10]-y[3]*h-y[5]*d*sqrt(2)/2;
	      auto v4 = y[10]-y[3]*h+y[5]*d*sqrt(2)/2;
  
	      auto w1 = y[11]-y[3]*d*sqrt(2)/2+y[4]*d*sqrt(2)/2;
	      auto w2 = y[11]-y[3]*d*sqrt(2)/2+y[4]*d*sqrt(2)/2;
	      auto w3 = y[11]+y[3]*d*sqrt(2)/2+y[4]*d*sqrt(2)/2;
	      auto w4 = y[11]+y[3]*d*sqrt(2)/2-y[4]*d*sqrt(2)/2;
              
              /* computation of the four aerodynamical forces in the body frame */
              /* F_i^a=Om_i*K*((u_i,v_i,w_i)-RT*Wa */
              /* calcul de RT*Wa d'abord */
              /* for now Wa = 0 */
              static const double Wa1 = 0;
              static const double Wa2 = 0;
              static const double Wa3 = 0;
              auto RTWax = RT11*Wa1+RT12*Wa2+RT13*Wa3;
              auto RTWay = RT21*Wa1+RT22*Wa2+RT23*Wa3;
              auto RTWaz = RT31*Wa1+RT32*Wa2+RT33*Wa3;
              
              auto F1x = Om1*(K11*(u1-RTWax)+K12*(v1-RTWay)+K13*(w1-RTWaz));
              auto F1y = Om1*(K21*(u1-RTWax)+K22*(v1-RTWay)+K23*(w1-RTWaz));
              auto F1z = Om1*(K31*(u1-RTWax)+K32*(v1-RTWay)+K33*(w1-RTWaz));
              
              auto F2x = Om2*(K11*(u2-RTWax)+K12*(v2-RTWay)+K13*(w2-RTWaz));
              auto F2y = Om2*(K21*(u2-RTWax)+K22*(v2-RTWay)+K23*(w2-RTWaz));
              auto F2z = Om2*(K31*(u2-RTWax)+K32*(v2-RTWay)+K33*(w2-RTWaz));
              
              auto F3x = Om3*(K11*(u3-RTWax)+K12*(v3-RTWay)+K13*(w3-RTWaz));
              auto F3y = Om3*(K21*(u3-RTWax)+K22*(v3-RTWay)+K23*(w3-RTWaz));
              auto F3z = Om3*(K31*(u3-RTWax)+K32*(v3-RTWay)+K33*(w3-RTWaz));
              
              auto F4x = Om4*(K11*(u4-RTWax)+K12*(v4-RTWay)+K13*(w4-RTWaz));
              auto F4y = Om4*(K21*(u4-RTWax)+K22*(v4-RTWay)+K23*(w4-RTWaz));
              auto F4z = Om4*(K31*(u4-RTWax)+K32*(v4-RTWay)+K33*(w4-RTWaz));
              
              /* computation of the three aerodynamical moments in the body frame */
	      auto Max = -d*sqrt(2)/2*F1z-h*F1y-d*sqrt(2)/2*F2z-h*F2y+d*sqrt(2)/2*F3z-h*F3y+d*sqrt(2)/2*F4z-h*F4y;
	      auto May = h*F1x-d*sqrt(2)/2*F1z+h*F2x+d*sqrt(2)/2*F2z+h*F3x+d*sqrt(2)/2*F3z+h*F4x-d*sqrt(2)/2*F4z;
	      auto Maz = d*sqrt(2)/2*F1y+d*sqrt(2)/2*F1x-d*sqrt(2)/2*F2y+d*sqrt(2)/2*F2x-d*sqrt(2)/2*F3y-d*sqrt(2)/2*F3x+d*sqrt(2)/2*F4y-d*sqrt(2)/2*F4x;

              //std:cout << getAAF(cmd_p).convert_int() << std::endl;
              
              /* controlled moments and vertical force are corrected using the aerodynamical effects */
              auto Mx = ((4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_r + (-4*C1*C1*Ct*d)*cmd_p*cmd_y)+Max;
              auto My = (-4*C1*C1*Ct*d*cmd_r*cmd_y + (4*Ct*d*thrust*C1*C1 + 4*C2*Ct*d*C1)*cmd_p)+May;
              auto Mz = (-2*C1*C1*Cd*cmd_r*cmd_p + (8*Cd*thrust*C1*C1 + 8*C2*Cd*C1)*cmd_y)+Maz;
              auto F  = Ct*C1*C1 *(cmd_p*cmd_p + cmd_r*cmd_r + 4.0*cmd_y*cmd_y + 4.0*thrust*thrust) + 8*Ct*C1*C2*thrust + 4*Ct*C2*C2+F1z+F2z+F3z+F4z;
              
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
              
              // derivatives of body speed u , v and w, corrected with the aerodynamical forces
              yp[9] = y[5]*y[10] - y[4]*y[11] + g*sinPitch+(F1x+F2x+F3x+F4x)/m;
              yp[10] = y[3]*y[11] - y[5]*y[9] - g*cosPitch*sinRoll+(F1y+F2y+F3y+F4y)/m;
              yp[11] = y[4]*y[9] - y[3]*y[10] + F/m - g*cosPitch*cosRoll;
              
              // //Z derivative coordinate
              // yp[12] = cosPitch*cosRoll*y[11] - sinPitch*y[9] + cosPitch*sinRoll*y[10];
              // // Z integrale for thrust setpoint calculation
              // yp[13] = velZ_sp - y[11];
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
        else if (syschoice == 6)  // self-driving car: parameters are not in the Jacobian
        {
            yp[0] = y[1];
            yp[1] = -params[0] *(y_prev[0] - 1.0) - params[1]*y_prev[1];   // pr = 1 is the reference position
        }
        else if (syschoice == 8)  // self-driving car with uncertain (but constant) coefficients
        {
            yp[0] = y[1];
            yp[1] = -param_inputs[0] *(y_prev[0] - 1.0) - param_inputs[1]*y_prev[1];   // pr = 1 is the reference position
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
            Jp[1][0] = - fullinputs[0]*J_prev[0][0]   - fullinputs[1]*J_prev[1][0];
            Jp[1][1] = - fullinputs[0]*J_prev[0][1]   - fullinputs[1]*J_prev[1][1];
            Jp[1][2] = -x_prev[0] - fullinputs[0]*J_prev[0][2] + 1 - fullinputs[1]*J_prev[1][2];
            Jp[1][3] = - fullinputs[0]*J_prev[0][3] - x_prev[1] - fullinputs[1]*J_prev[1][3];
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
