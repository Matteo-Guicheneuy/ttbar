// ************************************************************************* //
// Tools for std::string manipulations                                            //
//                                                                           //
// By Benjamin Fuks - 31.05.2012.                                            //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      // Mathematical functions         //
#include <math.h>
#include <string>     // Strings                        //
#include <sstream>    // String streams                 //
#include <fstream>    // File streams                   //
#include <cstring>    // In/Out streams                 //
#include <iostream>
#include <complex>
// -------- Classes ----------------------------------- //
#include "messages.h"        // Message services        //
#include "LHAPDF/LHAPDF.h"   // LHAPDF                  //

// ---------------------------------------------------- //

// ************************************************************************* //
// Printing a xsection                                                       //
// ************************************************************************* //
void DisplayXsec(const double &r, const double &e, const std::string &tag)
{
  std::ostringstream osr; osr << r; std::string rstr= osr.str();
  std::ostringstream ose; ose << e; std::string estr= ose.str();
  double prec=0.;
  if(r!=0.) prec=std::abs(100.*e/r);
  std::ostringstream opr; opr.precision(2); 
  opr << std::fixed << prec; std::string prstr= opr.str();
  info(tag + " results: " + rstr + " +/- " + estr + " pb  (@ " + prstr + "%)");
}


// ************************************************************************* //
// Other tools                                                               //
// ************************************************************************* //

double derivk(double (* f)(double x, int i0, const LHAPDF::PDF* F), double x, int i0, const LHAPDF::PDF* F ,double eps, int k)
{
    switch(k){
      case 1: return (3.*f(x-4.*eps,i0,F)-32.*f(x-3.*eps,i0,F)+168.*f(x-2.*eps,i0,F)-672.*f(x-eps,i0,F)+672.*f(x+eps,i0,F)-168.*f(x+2.*eps,i0,F)+32.*f(x+3.*eps,i0,F)-3.*f(x+4.*eps,i0,F))*pow(840.*eps,-1);
        break;
      case 2: return (-9.*f(x-4.*eps,i0,F)+128.*f(x-3.*eps,i0,F)-1008.*f(x-2.*eps,i0,F)+8064.*f(x-eps,i0,F)-14350*f(x,i0,F)+8064.*f(x+eps,i0,F)-1008.*f(x+2.*eps,i0,F)+128.*f(x+3.*eps,i0,F)-9.*f(x+4.*eps,i0,F))/5040.*pow(eps,-2);
        break;
      case 3: return (-7.*f(x-4.*eps,i0,F)+72.*f(x-3.*eps,i0,F)-338.*f(x-2.*eps,i0,F)+488.*f(x-eps,i0,F)-488.*f(x+eps,i0,F)+338.*f(x+2.*eps,i0,F)-72.*f(x+3.*eps,i0,F)+7.*f(x+4.*eps,i0,F))/240.*pow(eps,-3);
        break;
      default: return 0.;
        break;
    }
}

double derive_x_k(double (* f)(double x, int i0, const LHAPDF::PDF* F), double x, int i0, const LHAPDF::PDF* F, double epsilon,int k)
{ //k: order of derivative
  //For derivativ:
  //n=-1 backwar 
  //n=0 central 
  //n=1 forward
  double eps=epsilon;
  switch(k)
  {
    case 1: return f(x, i0, F)+x*derivk(f, x, i0, F,eps, 1);
      break;
    case 2: return f(x, i0, F)+3.*x*derivk(f, x, i0, F,eps, 1)+pow(x,2)*derivk(f, x, i0, F, eps, 2);
      break;
    case 3: return f(x, i0, F)+7.*x*derivk(f, x, i0, F, eps, 1)+6.*pow(x,2)*derivk(f, x, i0, F, eps, 2)+pow(x,3)*derivk(f, x, i0, F, eps, 3);
      break;
    default: return 0.;
      break;
  }

}
