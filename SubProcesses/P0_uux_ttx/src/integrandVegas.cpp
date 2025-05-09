// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
#include <fstream>
#include <cstddef>
#include <math.h>
#include <complex>
#include <boost/math/special_functions/zeta.hpp>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
// -------- Classes ----------------------------------- //
#include "process.h"	// Process                        //
#include "integrandDiff.h"
// -------- Functions --------------------------------- //
std::complex<double> N(double);
std::complex<double> Psi(std::complex<double>);
double Xcos(double);
double Xsin(double);

struct my_f_params { double sc; int mel; double M2; const LHAPDF::PDF* F;};  

extern "C"{
  //void update_mur_(double&);
  void ml5_0_sloopmatrix_thres_(double[4][4],complex<double>(*),double&,double[2][2],int&);
}

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern int chan;

//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//
double TotVegas(double *x, size_t dim, void *params)
{
  double res=0., bet, t, hgg,  u, Xbet;
  std::complex<double> temp=0, i(0.0,1.0), n, Nb;
  
  struct my_f_params * fparams = (struct my_f_params *)params;  
  double M2=fparams->M2, sc=fparams->sc; 
  const LHAPDF::PDF* F=fparams->F;
  int mel=fparams->mel;
  n=N(x[0]);
  Nb=n*exp(-Psi(1.));
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  //t=-M2/2.*(1-Xbet);

  if(mel==0)
    {
      std::complex<double> tmp=0.;      
      tmp=TraceBornDiff(n,Xbet,chan,M2)*2.;  res*=2./2./M_PI;
      res=std::imag(tmp*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
    }
    else if(mel==4)
    {
      std::complex<double> tmp=0.;
      double a=-1.;
      int k=2;
      tmp=TraceBornDiff(n+a,Xbet,chan,M2)*2.; // *2 because of qqb <-> qbq symmetry 
      std::complex<double> g_12=std::pow(n+a-1.,-2.*k);// fake pdf
      res=std::imag(tmp*g_12*xPDFDiff(x,chan,sc,M2,F,k)*GlobalDiffPDF(x,a,sc,M2));
        }

  else if(mel==1) //Diff
    {
      int Len=chan+2;
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));

      res-=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);

      tu=1;
      res+=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx2,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));

      res-=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx2,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
    }
    else if(mel==5) //diff pdf trick
    {
      res=0.;
    }

    else if(mel==2) //NLL 
    {
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      
      res+=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx2,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      
    }
    else if(mel==6)// NLL PDF trick
    {
      // Kinematics
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      // PDf trick variables
      double a=3.;
      int k=2;
      std::complex<double> g_12=std::pow(n+a-1.,-2.*k);// fake pdf

      ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);
      res=std::imag(TraceHSDiff(n+a,Xbet,chan,M2,xx,tu)*ColinearDiff(n+a,chan,M2)*g_12*xPDFDiff(x,chan,sc,M2,F,k)*GlobalDiffPDF(x,a,sc,M2));
      //std::cout << "res " <<TraceHSDiff(n+a,Xbet,chan,M2,xx,tu)*ColinearDiff(n+a,chan,M2) << std::endl;
      //std::cout  << g_12 << " "<< xPDFDiff(x,chan,sc,M2,F,k) << " "<< GlobalDiffPDF(x,a,sc,M2) << std::endl;
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      res+=std::imag(TraceHSDiff(n+a,Xbet,chan,M2,xx2,tu)*ColinearDiff(n+a,chan,M2)*g_12*xPDFDiff(x,chan,sc,M2,F,k)*GlobalDiffPDF(x,a,sc,M2));

    }

    else if(mel==3) //Exp
    {
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));

      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      
      res+=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx2,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
    }
    else if(mel==7) //Exp
    {
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;

      double a=3.;
      int k=2;
      std::complex<double> g_12=std::pow(n+a-1.,-2.*k);// fake pdf
      
      //ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag((TraceHSExpDiff(n+a,Xbet,chan,M2,xx,tu)+(1.+ColinearExpDiff(n+a,chan,M2))*TraceBornDiff(n+a,Xbet,chan,M2))*g_12*xPDFDiff(x,chan,sc,M2,F,k)*GlobalDiffPDF(x,a,sc,M2));

      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      
      res+=std::imag((TraceHSExpDiff(n+a,Xbet,chan,M2,xx2,tu)+(1.+ColinearExpDiff(n+a,chan,M2))*TraceBornDiff(n+a,Xbet,chan,M2))*g_12*xPDFDiff(x,chan,sc,M2,F,k)*GlobalDiffPDF(x,a,sc,M2));
    }

  
  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  return res;
}
