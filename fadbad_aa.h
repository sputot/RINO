/* ============================================================================
 File   : fadbad_aa.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This interface file allows us to use automatic differentiation from the FADBAD++ package on  intervals from FILIB++
 and affine forms AAF from the aaflib library
 ============================================================================ */


#ifndef FADBAD_AA_H
#define FADBAD_AA_H

#include "network_handler.h"

#include "filib_interval.h"

#include "aa_aaf.h"
#include "fadiff.h"
//#include "utils.h"

using namespace fadbad;

#define tol_noise 0.00001


template <> struct Op<interval>
{
    
    typedef interval Base;
    
    static Base myInteger(const int i) { return Base(double(i)); }
    static Base myZero() { return myInteger(0); }
    static Base myOne() { return myInteger(1);}
    static Base myTwo() { return myInteger(2); }
    
    static Base myPos(const Base& x) { return +x; }
    static Base myNeg(const Base& x) { return -x; }
    
    static Base& myCadd(Base& x, const Base& y) { return x+=y; }
    static Base& myCsub(Base& x, const Base& y) { return x-=y; }
    static Base& myCmul(Base& x, const Base& y) { return x*=y; }
    static Base& myCdiv(Base& x, const Base& y) { return x/=y; }
    static Base myInv(const Base& x) { return myOne()/x; }
    static Base mySqr(const Base& x) { return sqr(x); }
    
    static Base myPow(const Base& x, const int n) { return pow(x,n); }
    static Base myPow(const Base& x, const Base& y) { return pow(x,y); }
    
    static Base mySqrt(const Base& x) { return sqrt(x); }
    static Base myLog(const Base& x) { return log(x); }
    static Base myExp(const Base& x) { return exp(x); }
    static Base mySin(const Base& x) { return sin(x); }
    static Base myCos(const Base& x) { return cos(x); }
    static Base myTan(const Base& x) { return tan(x); }
    static Base myAsin(const Base& x) { return asin(x); }
    static Base myAcos(const Base& x) { return acos(x); }
    static Base myAtan(const Base& x) { return atan(x); }
    
    static Base mySig(const Base& x) { /*cout << "est bien appele (interv)";*/ return  1./(1.+exp(-x));} // act_sigmoid(x); } return act_sigmoid(x); }//
    static Base myTanh(const Base& x) { cout << "est bien appele aussi (interv)"; return  2.0/(1.+ exp(-2.0*x)) - 1.;} // ::act_tanh(x); }
    
    static bool myEq(const Base& x, const Base& y) { return x==y; }
    static bool myNe(const Base& x, const Base& y) { return x!=y; }
};


template <> struct Op<AAF>
{
  
  typedef AAF Base;

  static Base myInteger(const int i) { return Base(double(i)); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }

  static Base myPos(const Base& x) { return  x; }
  static Base myNeg(const Base& x) { return  -x;  }
  
  static Base& myCadd(Base& x, const Base& y) { return x+=y;   }
  static Base& myCsub(Base& x, const Base& y) { return x-=y;   }
  static Base& myCmul(Base& x, const Base& y) { return x*=y;  }
  static Base& myCdiv(Base& x, const Base& y) { return x=x/y; }
  static Base myInv(const Base& x) { return myOne()/x;  }

  static Base mySqr(const Base& x) { return sqr(x); }
	
  static Base myPow(const Base& x, const int n) { return x^n; }
  static Base myPow(const Base& x, const Base& y) { return x^y; }

  static Base mySqrt(const Base& x) { return sqrt(x);  }
  static Base myLog(const Base& x) { return log(x);  }
  static Base myExp(const Base& x) { return exp(x);  }
  static Base mySin(const Base& x) { return sin(x);  }
  static Base myCos(const Base& x) { return cos(x);  }
  static Base myTan(const Base& x) { return tan(x);  }
  static Base myAsin(const Base& x) { return myInteger(0); throw std::invalid_argument( "asin non prevu");}
  static Base myAcos(const Base& x) { return myInteger(0); throw std::invalid_argument( "acos non prevu"); }
  static Base myAtan(const Base& x) { return atan(x); }
    
    static Base mySig(const Base& x) { /*cout << "est bien appele"; */AAF y = x; y.sumup(tol_noise); return act_sigmoid(y); }// 1./(1.+exp(-x));} //act_sigmoid(x); }
    static Base myTanh(const Base& x) { cout << "est bien appele aussi"; AAF y = x; y.sumup(tol_noise); return act_tanh(y); return 2.0/(1.+ exp(-2.0*x)) - 1.;} // ::act_tanh(x); }
  
  /*  static Base myPos(const Base& x) { return +x; }
  static bool myNe(const Base& x, const Base& y) { return x!=y; }
  */ 
  static bool myEq(const Base& x, const Base& y) { return x==y; }
  static bool myNe(const Base& x, const Base& y) { return !(x==y); }
};

#endif
