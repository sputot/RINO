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



void define_system_dim(int argc, char* argv[]);  // define the dimensions of your system (ODE or DDE)



// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &params_inputs, vector<AAF> &param_inputs_center, vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for DDEs : functions that initialize the DDE on first time period
vector<T<AAF>> Initfunc(const T<AAF>& t, vector<AAF> &beta);
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta);

// defining analytical solution if any for comparison
vector <interval> AnalyticalSol(double t, vector<AAF> &beta, double d0);

// reading sysdim, jacdim, etc
void readfromfile_system_dim(const char * params_filename, int &sysdim, int &jacdim, int &sysdim_params, int &nb_subdiv_init);

// d0 and t_begin and nb_subdiv are for DDEs only, rest are common to ODE and DDE
void read_parameters(const char * params_filename, double &tau, double &t_end, double &d0, double &t_begin, int &order, int &nb_subdiv);

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(int argc, char* argv[], double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide);



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
          else if (syschoice == 8)
          {
              yp[0] = - y[0]*y[0]*y[0];
          }
          else if (syschoice == 9)  // self-driving car with constant parameters
          {
              yp[0] = y[1];
              yp[1] = -param_inputs[0] *(y[0] - 1.0) - param_inputs[1]*y[1];   // pr = 1 is the reference position
          }
          else if (syschoice == 10)  // self-driving car with constant parameters
          {
              yp[0] = y[1];
              yp[1] = -y[2] *(y[0] - 1.0) - y[3]*y[1];   // pr = 1 is the reference position
              yp[2] = 0; // constant parameter Kp
              yp[3] = 0; // constant parameter Kd
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
