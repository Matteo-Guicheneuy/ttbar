// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
#include <cstddef>
#include <math.h>
#include <complex>
#include <cstdlib>

// -------- Classes ----------------------------------- //

#include "process.h"	// Process//
#include "LHAPDF/LHAPDF.h"

extern "C" {
  double ct18pdf_(int&, double&, double&);              //            
  double ctq6pdf_(int&, double&, double&);              //            
}; 

// ---------------------------------------------------- //
extern double MZ, muF, muR, alpha, tw, gu, gd, alpha_s, beta0;
extern int Nf;
using namespace std;
using namespace LHAPDF;

complex<double> Psi(complex<double> N)
{
  complex<double> res(0.,0.);
  while(real(N)<10.) { res=res-1./N; N=N+1.; }
  res=res+log(N)-1./(2.*N)-pow(N,-2.)/2520.*(210.+pow(N,-2.)*(-21.+pow(N,-2.)*10.));
  return res;
}

void SetLHAPDF(const PDF* F, double &x, double q[5], double qb[5], double &g)
{

  for(int i0=-5; i0<6; i0++){
    if(i0<0) qb[-1-i0]=F->xfxQ2(i0,x,muF*muF)/x;
    else if(i0==0) g=F->xfxQ2(i0,x,muF*muF)/x;
    else q[i0-1]=F->xfxQ2(i0,x,muF*muF)/x;
  }
  
}

double FittingPDF(double &x, double A[8][8], int ID)
{
  double res=0;
  res=A[ID][0]*pow(x,A[ID][1])*pow(1.-x,A[ID][2])*(1.+A[ID][3]*pow(x,0.5)+x*A[ID][4]+pow(x,1.5)*A[ID][5]+pow(x,2)*A[ID][6]+pow(x,2.5)*A[ID][7]);
  return res;
}

void SetPDFfit(double &x, double A[8][8], double q[5], double qb[5], double &g)
{
  int ID=0;
  g=FittingPDF(x,A,ID);
  ID=1;
  q[0]=FittingPDF(x,A,ID);
  ID=2;
  q[1]=FittingPDF(x,A,ID);
  ID=3;
  qb[0]=FittingPDF(x,A,ID);
  ID=6;
  qb[1]=FittingPDF(x,A,ID);
  ID=4;
  q[2]=FittingPDF(x,A,ID);
  qb[2]=FittingPDF(x,A,ID);
  ID=5;
  q[4]=FittingPDF(x,A,ID);
  qb[4]=FittingPDF(x,A,ID);
  ID=7;
  q[3]=FittingPDF(x,A,ID);
  qb[3]=FittingPDF(x,A,ID);
}

complex<double> Gamma(const complex<double> x) {
    const int intx = (int)real(x);
    const double n = -(double)(intx < 0 ? -intx : intx);
    if (real(x) == n && imag(x) == 0.0) {
        cout << "Gamma(" << x << ") undefined\n";
        exit(0);
    }

    // Works with Re(xx) > 0
    const complex<double> xx = (real(x) > 0.0 ? x : -x);

    // Magic numbers for Gamma function
    const double q0 = 75122.6331530;
    const double q1 = 80916.6278952;
    const double q2 = 36308.2951477;
    const double q3 = 8687.24529705;
    const double q4 = 1168.92649479;
    const double q5 = 83.8676043424;
    const double q6 = 2.50662827511;

    complex<double> gamma = exp((xx + .5) * log(xx + 5.5) - xx - 5.5) *
                            (q0 + q1 * xx + q2 * pow(xx, 2) + q3 * pow(xx, 3) + q4 * pow(xx, 4) +
                             q5 * pow(xx, 5) + q6 * pow(xx, 6))
                            / xx / (xx + 1.0) / (xx + 2.0) / (xx + 3.0) / (xx + 4.0) / (xx + 5.0) / (xx + 6.0);

    return (x == xx ? gamma : -M_PI / xx / sin(M_PI * xx) / gamma);
}

complex<double> B2(complex<double> a, double b){ // Gamma(a)Gamma(b)/Gamma(a+b)
  complex<double> res=0;
  res=Gamma(a)*Gamma(b)/Gamma(a+b);
  
  return res;
}


complex<double> Fbis(complex<double> &N, double A[8][8], int &i0)
{
  complex<double> res;
  res=B2(A[i0][1]+N,A[i0][2]+1.)+A[i0][3]*B2(A[i0][1]+N+0.5,A[i0][2]+1.)+A[i0][4]*B2(A[i0][1]+N+1.,A[i0][2]+1.)+A[i0][5]*B2(A[i0][1]+N+1.5,A[i0][2]+1.)+A[i0][6]*B2(A[i0][1]+N+2.,A[i0][2]+1.)+A[i0][7]*B2(A[i0][1]+N+2.5,A[i0][2]+1.);

  res*=A[i0][0];
  
  return res;
}

complex<double> B(complex<double> &n, double a, double b){ // Gamma(N+a)/Gamma(N+b+a)
  complex<double> res=0;


  res=(1. - (b*(-1.+ + 2*a + b))/(2.*n) + (b*(1. + b)*(2. + 12*pow(a,2) + 12*a*(-1. + b) - 5*b + 3*pow(b,2)))/(24.*pow(n,2)) - (b*(2. + 3*b + pow(b,2))*(8*pow(a,3) + 12*pow(a,2)*(-1. + b) + pow(-1. + b,2)*b + 2*a*(2. - 5*b + 3*pow(b,2))))/(48.*pow(n,3)) + (b*(6. + 11*b + 6*pow(b,2) + pow(b,3))*(-8. + 240*pow(a,4) + 480*pow(a,3)*(-1. + b) + 18*b + 120*a*pow(-1. + b,2)*b + 5*pow(b,2) -  30*pow(b,3) + 15*pow(b,4) + 120*pow(a,2)*(2. - 5*b + 3*pow(b,2))))/(5760.*pow(n,4)) -   (b*(24. + 50*b + 35*pow(b,2) + 10*pow(b,3) + pow(b,4))* (96*pow(a,5) + 240*pow(a,4)*(-1. + b) + 120*pow(a,2)*pow(-1. + b,2)*b + 80*pow(a,3)*(2. - 5*b + 3*pow(b,2)) +      pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 2*a*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/(11520.*pow(n,5)) +      (b*(120. + 274*b + 225*pow(b,2) + 85*pow(b,3) + 15*pow(b,4) + pow(b,5))*        (96. + 4032*pow(a,6) + 12096*pow(a,5)*(-1. + b) - 236*b + 10080*pow(a,3)*pow(-1. + b,2)*b - 84*pow(b,2) + 539*pow(b,3) -           315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6) + 5040*pow(a,4)*(2. - 5*b + 3*pow(b,2)) +           252*a*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 252*pow(a,2)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/   (2.90304e6*pow(n,6)) - (b*(720. + 1764*b + 1624*pow(b,2) + 735*pow(b,3) + 175*pow(b,4) + 21*pow(b,5) + pow(b,6))*
        (1152*pow(a,7) + 4032*pow(a,6)*(-1. + b) + 5040*pow(a,4)*pow(-1. + b,2)*b + 2016*pow(a,5)*(2. - 5*b + 3*pow(b,2)) +    252*pow(a,2)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +           pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +       168*pow(a,3)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +        2*a*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/(5.80608e6*pow(n,7)) +    (b*(5040. + 13068*b + 13132*pow(b,2) + 6769*pow(b,3) + 1960*pow(b,4) + 322*pow(b,5) + 28*pow(b,6) + pow(b,7))*    (-1152. + 34560*pow(a,8) + 138240*pow(a,7)*(-1. + b) + 3088*b + 241920*pow(a,5)*pow(-1. + b,2)*b + 884*pow(b,2) -      8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8) +       80640*pow(a,6)*(2. - 5*b + 3*pow(b,2)) + 20160*pow(a,3)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +        240*a*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         10080*pow(a,4)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +      240*pow(a,2)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/   (1.3934592e9*pow(n,8)) - (b*(40320. + 109584*b + 118124*pow(b,2) + 67284*pow(b,3) + 22449*pow(b,4) + 4536*pow(b,5) +        546*pow(b,6) + 36*pow(b,7) + pow(b,8))*(7680*pow(a,9) + 34560*pow(a,8)*(-1. + b) + 80640*pow(a,6)*pow(-1. + b,2)*b +        23040*pow(a,7)*(2. - 5*b + 3*pow(b,2)) + 10080*pow(a,4)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +         240*pow(a,2)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         4032*pow(a,5)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +       pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +         160*pow(a,3)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +         2*a*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) +            135*pow(b,8))))/(2.7869184e9*pow(n,9)) + (b*(362880. + 1026576*b + 1172700*pow(b,2) + 723680*pow(b,3) + 269325*pow(b,4) + 63273*pow(b,5) + 9450*pow(b,6) +         870*pow(b,7) + 45*pow(b,8) + pow(b,9))*(7680. + 101376*pow(a,10) + 506880*pow(a,9)*(-1. + b) - 22112*b +      1520640*pow(a,7)*pow(-1. + b,2)*b - 3960*pow(b,2) + 62524*pow(b,3) - 56958*pow(b,4) - 1265*pow(b,5) + 20559*pow(b,6) -       5082*pow(b,7) - 1980*pow(b,8) + 495*pow(b,9) + 99*pow(b,10) + 380160*pow(a,8)*(2. - 5*b + 3*pow(b,2)) +         266112*pow(a,5)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +          10560*pow(a,3)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +        88704*pow(a,6)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +           132*a*pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +           5280*pow(a,4)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +   132*pow(a,2)*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8))))/(3.678732288e11*pow(n,10)))/pow(n,b);

  
  if (res != res) std::cout << "Nan in B" << std::endl;

  if (abs(res) > 1e100) std::cout << "B = " << res << " a= " << a << " b= " << b << std::endl;
  return res;
}


complex<double> F3(complex<double> &N, double A[8][8], int &i0) // Calcul simplifié dans le cas |N| > 1e2
{
  complex<double> res=0;

  res= B(N,A[i0][1],1.+A[i0][2])+B(N,A[i0][1]+0.5,1.+A[i0][2])*A[i0][3]+B(N,A[i0][1]+1.,1.+A[i0][2])*A[i0][4]+B(N,A[i0][1]+1.5,1.+A[i0][2])*A[i0][5]+B(N,A[i0][1]+2.,1.+A[i0][2])*A[i0][6]+B(N,A[i0][1]+2.5,1.+A[i0][2])*A[i0][7];
  res= res*A[i0][0]*Gamma(A[i0][2]+1.);
  
  return res;
}


void SetPDFN(complex<double> &N, complex<double> q[5], complex<double> qbar[5], complex<double> &g, double A[8][8])
{

  if(abs(N) < 140.){
    
    int i0=0; g=Fbis(N,A,i0);//gluon
    
    i0=5; qbar[4]=Fbis(N,A,i0); //bbar -> b from sea
    q[4]=Fbis(N,A,i0); //b -> b from sea
    
    i0=7; qbar[3]=Fbis(N,A,i0); //cbar -> c from sea
    q[3]=Fbis(N,A,i0); //c -> c from sea
  
    i0=4; qbar[2]=Fbis(N,A,i0); //sbar -> s from sea
    q[2]=Fbis(N,A,i0); //s -> s from sea
    
    i0=6; qbar[1]=Fbis(N,A,i0); //ubar -> u from sea
    i0=2; q[1]=Fbis(N,A,i0); //u -> u from valence  +F(N,A,6) taken into account in fit
  
    i0=3; qbar[0]=Fbis(N,A,i0); //dbar -> d from sea
    i0=1; q[0]=Fbis(N,A,i0); //d -> d from valence  +Fbis(N,A,3) taken into account in fit
  }

  else
    {
      int i0=0; g=F3(N,A,i0);//gluon
      //std::cout << "N=" << N << ", F(N)=" << g << std::endl;
      
      i0=5; qbar[4]=F3(N,A,i0); //bbar -> b from sea
      q[4]=F3(N,A,i0); //b -> b from sea
      
      i0=7; qbar[3]=F3(N,A,i0); //cbar -> c from sea
      q[3]=F3(N,A,i0); //c -> c from sea
      
      i0=4; qbar[2]=F3(N,A,i0); //sbar -> s from sea
      q[2]=F3(N,A,i0); //s -> s from sea
      
      i0=6; qbar[1]=F3(N,A,i0); //ubar -> u from sea
      i0=2; q[1]=F3(N,A,i0) ; //u -> u from valence  +F(N,A,6) taken into account in fit
      
      i0=3; qbar[0]=F3(N,A,i0); //dbar -> d from sea
      i0=1; q[0]=F3(N,A,i0); //d -> d from valence  +Fbis(N,A,3) taken into account in fit
      
    }
 

}


// complex<double> EvolOperator(complex<double> &N, complex<double> &EigenV, double Q)
// {
//   complex<double> res(0.0,0.0), Nb=N*exp(-Psi(1.));
  
//   res=pow((1.+beta0/M_PI*alpha_s*log(Q/Nb/muR))/(1.+beta0/M_PI*alpha_s*log(muF/muR)),EigenV/beta0);
//   return res;
// }

// void EvolvePDF(complex<double> &N, complex<double> q[5], complex<double> qbar[5], complex<double> &g,double M2)
// {
//   double nf=(double)Nf, cmp=0;
//   complex<double> V[5], NS[5], Sigma(0.,0.), tmp(0.0,0.0), Gqq, Ggg, Gqg, Ggq, rp, rm;
//   complex<double> VE[5], NSE[5], SigmaE(0.,0.);
//   double Q=pow(M2,0.5);
    
//   Gqq=4./3.*(3./2.+1./(N*(N+1.))-2.*(Psi(N+1.)-Psi(1.)));
//   Gqg=(2.+N+N*N)/(2.*N*(N+1.)*(N+2.));
//   Ggq=4./3.*(2.+N+N*N)/(N*(N*N-1.));
//   Ggg=beta0+2.*3.*(1./(N*(N-1.))+1./((N+1.)*(N+2.))+Psi(1.)-Psi(N+1.));
		   
//   rp=0.5*(Ggg+Gqq+pow(pow(Gqq-Ggg,2.)+8.*nf*Gqg*Ggq,0.5));
//   rm=0.5*(Ggg+Gqq-pow(pow(Gqq-Ggg,2.)+8.*nf*Gqg*Ggq,0.5));
  
  
//   for(int i0=0; i0<Nf; i0++){
//     cmp+=1.;
    
//     V[i0]=q[i0]-qbar[i0];
//     Sigma+=q[i0]+qbar[i0];
//     NS[i0]=Sigma-cmp*(q[i0]+qbar[i0]);
    
//     VE[i0]=V[i0]*EvolOperator(N,Gqq,Q);
//     NSE[i0]=NS[i0]*EvolOperator(N,Gqq,Q);
//   }

//   SigmaE=EvolOperator(N,rp,Q)*((Gqq-rm)*Sigma+2.*nf*Gqg*g)/(rp-rm)-EvolOperator(N,rm,Q)*((Gqq-rp)*Sigma+g*Gqg*2.*nf)/(rp-rm);

//   //Evolved gluon PDF
//   g=EvolOperator(N,rp,Q)*((Ggg-rm)*g+Ggq*Sigma)/(rp-rm)-EvolOperator(N,rm,Q)*((Ggg-rp)*g+Ggq*Sigma)/(rp-rm);

//   //Evolved quark PDF
//   for (int i0 = Nf-1; i0 >= 0; i0--) {  
//      q[i0] = (1.0 / nf * SigmaE - 1.0 / cmp * NSE[i0] + tmp + VE[i0]) * 0.5;
//      qbar[i0] = (1.0 / nf * SigmaE - 1.0 / cmp * NSE[i0] + tmp - VE[i0]) * 0.5;
//      tmp += 1.0 / cmp / (cmp - 1.0) * NSE[i0];
//      cmp -= 1.0;
//    }
  
// }
