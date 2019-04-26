/* ============================================================================
 File   : filib_interval.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This  file is an interface to the intervals of FILIB++ library
 ============================================================================ */


#ifndef FILIB_INTERVAL_H
#define FILIB_INTERVAL_H


#include <interval/interval.hpp> 
#include <iostream> 
#include <sstream> 
#include <string> 


   typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> interval;
   //   typedef filib::interval<double> interval;


inline double inf(const interval&a){return a.inf();}

inline double sup(const interval&a){
return a.sup();}

inline double midpoint(const interval&a){
return a.mid();
}

inline double width(const interval&a)
{return a.diam();}


inline double mag(const interval&a){return a.mag();}


inline bool subseteq(const interval&a,const interval&b)
{return filib::subset(a,b);}

inline bool interior(const interval&a,const interval&b)
{return filib::interior(a,b);}


inline bool disjoint(const interval&a,const interval&b)
{
return filib::disjoint(a,b);
}

inline bool intersect(interval&c,
const interval&a,const interval&b)
{
if(filib::disjoint(a,b))
return false;
c= filib::intersect(a,b);
return true;
}

inline interval hull(const interval&a,const interval&b)
{
    return filib::hull(a,b);
}

inline interval pi(){
  return interval::PI(); //filib::interval<double> ::PI();
}

inline interval empty(){
  return interval::EMPTY(); // filib::interval<double> ::EMPTY();
}

inline interval entire() {
return interval::ENTIRE(); 
}

inline interval pow(const interval&a,const interval&b){
return filib::pow(a,b);
}

inline interval pow(const interval&a,int b){
return filib::power(a,b);
}


inline interval exp(const interval&a){
return filib::exp(a);
}

inline interval log(const interval&a){
return filib::log(a);
}

inline interval sqr(const interval&a){
return filib::sqr(a);
}

inline interval sqrt(const interval&a){
return filib::sqrt(a);
}


inline interval sin(const interval&a){
return filib::sin(a);
}

inline interval cos(const interval&a){
return filib::cos(a);
}

inline interval tan(const interval&a){
return filib::tan(a);
}


inline interval asin(const interval&a){
return filib::asin(a);
}

inline interval acos(const interval&a){
return filib::acos(a);
}

inline interval atan(const interval&a){
return filib::atan(a);
}



inline interval string_to_interval(const char*s)
    {
        return interval(s,s);
    }




#endif
