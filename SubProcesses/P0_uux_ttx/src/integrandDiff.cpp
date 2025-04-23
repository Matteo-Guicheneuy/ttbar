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
#include "process.h"	// Process//
#include "setting.h"
#include "main.h"
#include "SoftFunctions.h"
// -------- Functions --------------------------------- //
double Sqrt(double);
std::complex<double> Psi(std::complex<double>);
std::complex<double> Gamma(const std::complex<double>);
std::complex<double> B(std::complex<double>&,double,double);
double Ci(int);
double B1(int);
std::complex<double> N(double);
double Xcos(double);
double Xsin(double);
std::complex<double> GammGlobal(double*, int, int, int, double, int);
double H0(std::complex<double>, double ,int,int,int,double);
void pdfFit(double&, double (*)[8], double, double, double);

extern "C"{
  //void update_mur_(double&);
  void ml5_0_sloopmatrix_thres_(double[4][4],complex<double>(*),double&,double[2][2],int&);
}

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern double A[8][8], A1min;
extern int Nc, Nf, chan;
extern string namePDF;


std::complex<double> TraceBornDiff(std::complex<double > N, double Xbet,int chan, double M2)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  std::complex<double> res=0.;    
  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=S0(ka,ji,chan)*H0(N,Xbet,ji,ka,chan,M2);
	 }
     }  
  return res;
}

std::complex<double> TraceHSDiff(std::complex<double > N,double Xbet ,int chan, double M2, complex<double> *xx, int tu)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  std::complex<double> res=1.;
  std::complex<double> **R, **Rm, **Rdag, **Rdagm1, *Lamb;
  std::complex<double> H0t[Len][Len]={0}, i(0.0,1.0);

  std::complex<double> Nb, DeltaProc;
  std::complex<double> St[Len][Len]={0};
  std::complex<double> S0t[Len][Len]={0};
  std::complex<double> H1t[Len][Len]={0};
  // Def of matrices for the diagonalisation of the Sudakov anomalous dimension matrix
  R=new std::complex<double> *[Len];
  Rm=new std::complex<double> *[Len];
  Rdag=new std::complex<double> *[Len];
  Rdagm1=new std::complex<double> *[Len];
  Lamb=new std::complex<double> [Len];// valeur propre
  
  for(int h=0; h < Len; h++)
    {
      R[h]=new std::complex<double>[Len];
      Rdag[h]=new std::complex<double>[Len];
      Rm[h]=new std::complex<double>[Len];
      Rdagm1[h]=new std::complex<double>[Len];
    }
  
  EigenvGlobal(N,Xbet,M2,chan,Rdag,Rdagm1,R,Rm,Lamb,tu); // Computes rotation matrices and eigenvalues from soft anomalous dimension matrix
  Nb=N*exp(-Psi(1.));
  
  double alpha_p=0.118;
  
  for(int k=0; k < Len; k++){
    for(int l=0; l < Len; l++){
      // Hard and Soft functions
	    for(int mm=0; mm < Len; mm++){
          for(int nn=0; nn < Len; nn++){
            // Firt Order Hard function
		        H0t[k][l]+=R[k][mm]*H0(N,Xbet,mm,nn,chan,M2)*Rdag[nn][l];
            // Second Order Hard function
		        H1t[k][l]+=R[k][mm]*xx[1+nn*4+mm*Len*4]*pow(alpha_s/alpha_p,3)/2./M2*Rdag[nn][l];

            // Firt Order Soft function
		        S0t[k][l]+=Rdagm1[k][mm]*S0(mm,nn,chan)*Rm[nn][l];
            // First + Second Order Soft function
		        St[k][l]+=Rdagm1[k][mm]*(S0(mm,nn,chan)+S1t(N,Xbet,mm,nn,chan,M2,tu))*Rm[nn][l];
          }
	      }
      // Exp of Soft anomalous dimension matrix: exp(ln(1-2lambda)/beta0*Gamma_S):
	    //DeltaProc=exp(-(Lamb[l]+std::conj(Lamb[k]))/2./beta0*log(1.-alpha_s/M_PI*beta0*log(Nb))); //HS convention
	    DeltaProc=exp((Lamb[l]+std::conj(Lamb[k]))/2./beta0*log(1.-alpha_s/M_PI*beta0*log(Nb))); //Me convention
      
      //Product with the Soft function
	    St[k][l]*=DeltaProc;
	    S0t[k][l]*=DeltaProc;
	  }
  }
  
  res=0.; 
  // Compute the trace of Hard and Soft
   for(int ka=0; ka < Len; ka++){
       for(int ji=0; ji < Len; ji++){
        // Porduct with the Hard function and Trace
	      res+=St[ka][ji]*H0t[ji][ka];
	      res+=S0t[ka][ji]*H1t[ji][ka];// the missing term at second order 
	      }
      }
    // kill pointeur
   for (int i=0; i<Len; i++)
     {
       delete [] R[i];
       delete [] Rm[i];
       delete [] Rdagm1[i];
       delete [] Rdag[i];
     }
    
   delete [] R;
   delete [] Rm;
   delete [] Rdag;
   delete [] Rdagm1;
   delete [] Lamb;
   
  
   if(isfinite(abs(res))==0){
     std::cout << "!! Result non finite in Res !!" << std::endl;
     exit(1);
   }
   
   return res;
 
}

std::complex<double> TraceHSExpDiff(std::complex<double > N, double Xbet ,int chan, double M2, complex<double> *xx, int tu)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  double alpha_p=0.118;
  
  std::complex<double> res=1.;

  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++){
       for(int ji=0; ji < Len; ji++){
	        res+=(S1t(N,Xbet,ka,ji,chan,M2,tu)+SG(N,Xbet,ka,ji,chan,M2,tu))*H0(N,Xbet,ji,ka,chan,M2);
	        res+=S0(ka,ji,chan)*xx[1+ka*4+ji*Len*4]*pow(alpha_s/alpha_p,3)/2./M2; //Virtual contribution
	      }
    }

  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in Tr expanded !!" << std::endl;
    exit(1);
  }
  
   return res;
}

// *************************************************** //
// Soft colinear radiation(x[0])s //
// *************************************************** //

std::complex<double> ColinearDiff(std::complex<double > N, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  //Change of variable
  Nb=N*exp(-Psi(1.));

  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  // g1, g2 and corresponding G functions of the N variable, independant of the correction considered (global) cf. Nucl. Phys. B. 529 (1998) Bociani, Catani, Mangano, Nason or 1005.2909 Debove, Fuks, Klasen or 0805.1885 
  // See HuaSheng note, thes are the non colinear imporve version
  g1=Ci(chan)/beta0/lambdN*(2.*lambdN+(1.-2.*lambdN)*log(1.-2.*lambdN));
  
  g2=Ci(chan)/beta0*log(1.-2.*lambdN)*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;
  g2+=-(2.*lambdN+log(1.-2.*lambdN))/pow(beta0,2)*Ci(chan)*(67.*double(Nc)/18.-double(Nc)/6.*pow(M_PI,2)-5.*double(Nf)/9.);
  g2+=beta1/pow(beta0,3)*Ci(chan)*(2.*lambdN+log(1.-2.*lambdN)+0.5*pow(log(1.-2.*lambdN),2));

  //g2+=g1*(-Psi(1.)); // N resummation
  //g2+=Psi(1.)/beta0*Ci(chan)/lambdN*(2.*lambdN+log(1.-2.*lambdN));

  //Sudakov Colinear= Delta_i^2 for qq > ttbar or gg > ttbar
  g1*=2.;
  g2*=2.;

  G=std::exp(g1*log(Nb)+g2);
  
  // result
  return G;
}

std::complex<double> ColinearExpDiff(std::complex<double > N, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  
  // Change of variable
  Nb=N*exp(-Psi(1.));
  
  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  // def of the sudakow factor in which we only keep the first order
  g1=2.*Ci(chan)*lambdN/beta0;
  
  g2=-2.*Ci(chan)/beta0*lambdN*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;

  //g2+=g1*(-Psi(1.)); // N resummation
  //g2+=(-Psi(1.))*lambdN*2.*Ci(chan)/beta0;
  
  g1*=2.;
  g2*=2.;
  // G expension
  G=1.+g1*log(Nb)+g2;
  
  if(isfinite(abs(G))==0){
    std::cout << "!! Result non finite in Colinear expanded !!" << std::endl;
    exit(1);
  }
  // result
  return G-1.;
}

std::complex<double> MellinPDFDiff(std::complex<double > N, int chan, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform           
{
  std::complex<double> i(0.0,1.0), fAB=0, res=0;

  std::complex<double> q[5], g, qbar[5], Npdf;
  Npdf=N;
  SetPDFN(Npdf,q,qbar,g,A); // PDFs with N


  // Symmetrized PDF depend on channel, 0 is qqbar and 1 is gg                                     
  if(chan==0)
    {
      //fAB=2.*Nf/(Npdf-1.)/(Npdf*Npdf)/(Npdf-1.); // Ansatz qqb    
      for(int j=0; j < Nf; j++) fAB+=q[j]*qbar[j];  // Carefull, beyond LO not symmetric in cos(theta) -> no factor 2 but 2 terms
      if(remove_charm_quark) fAB-=q[3]*qbar[3];
    }
  else if(chan==1) exit(1);       

  else std::cout << "Error in chan ID" << std::endl;

  //Global factors                                              
  res=fAB;
  return res;
}



std::complex<double> GlobalDiff(double *x, double sc, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform
{  
  std::complex<double> i(0.0,1.0), res=0, jac=0;
  double jac2=0, bet=0;
  double tau = M2/sc;
  
  // Change of variable
  jac = (cos(Phi)+i*sin(Phi))/x[0];
  
  // Change of variable 2
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  jac2=20.-10.*pow(tanh(20.*x[1]-5.),2)-10.*pow(tanh(20.*x[1]-15.),2);
  
  //Global factors
  res=pow(tau,-N(x[0])+1.)*2./2./M_PI; //*2 from Im part, /2*M_PI from inverse Mellin transform, N-1 because the functions are supposed to be evaluated at N+1 when taking the Mellin inverse in N
  res*=bet/16./M_PI/M2; 
  
  // Result
  res=res*jac*jac2*0.38937966e9;//*jac3;

  return res;
}
// *************************************************** //
// PDF Trick //
// *************************************************** // 
double xPDFDiff(double *x , int chan, double sc, double M2,const LHAPDF::PDF* F ,int k) // Resummed formula, the numericall integration does the inversed Mellin transform           
{
  double fAB=0., tau=M2/sc, a=-1.;
  double xa=tau/pow(tau,x[2]), xb=tau/pow(tau/xa,x[3])/xa;
  //PDF* F=mkPDF(namePDF,0);
  fAB=1./2.*Luminosity(xa,xb,F,a,usePDFAnzat,true,k); 
  return fAB;
}


std::complex<double> GlobalDiffPDF(double *x, double a,double sc, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform
{  
  std::complex<double> i(0.0,1.0), res=0, jac=0;
  double jac2=0., bet=0., jacpdf=0.;
  double tau = M2/sc;
  
  // Change of variable
  jac = (cos(Phi)+i*sin(Phi))/x[0];
  
  // Change of variable 2
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  jac2=20.-10.*pow(tanh(20.*x[1]-5.),2)-10.*pow(tanh(20.*x[1]-15.),2);
  
  //Global factors
  res=bet/16./M_PI/M2; 

  // PDF trick factor
  double xa=tau/pow(tau,x[2]), xb=tau/pow(tau/xa,x[3])/xa;
  jacpdf= xa*xb*pow(log(tau),2.)*x[2];
    //std::cout << x[2] <<" xa= " << xa << " xb= " << xb << " tau= " << tau <<" jacpdf = " << jacpdf << std::endl;
  // Result
  res*=pow(tau/xa/xb,-N(x[0])+1.-a)/xa/xb*2./2./M_PI; //*2 from Im part, /2*M_PI from inverse Mellin transform, N-1 because the functions are supposed to be evaluated at N+1 when taking the Mellin inverse in N
  res=res*jac*jac2*jacpdf*0.38937966e9;//*jac3;

  return res;
}
//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//
double TotDiff(double *x, double& sc, int& mel, double& M2,const LHAPDF::PDF* F)
{
  double res=0., bet, t, hgg,  u, Xbet;
  std::complex<double> temp=0, i(0.0,1.0), Nb;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  std::complex<double> n=N(x[0]);

  if(mel==0)
    {
      std::complex<double> tmp=0.;
      tmp=TraceBornDiff(n,Xbet,chan,M2)*2.; // *2 because of qqb <-> qbq symmetry 
      res+=1./M_PI*std::imag(tmp*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
    }

  else if(mel==1) //Diff and NLL
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
      // quadri implusion [nparticle][t,x,y,z]
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
      //xx[ 0 + nn*4+mm*4*Len ] est l'élément (nn,mm) de la matrice H0
      //xx[ 1 + nn*4+mm*4*Len ] est l'élément (nn,mm) de la matrice H1
      //xx[ 2 + nn*4+mm*4*Len ] est l'élément (nn,mm) de la matrice V-2 (poles en 1/eps^2 qui vient du calcul après renormalisation -> divergence IR à compenser avec les PDF)
      //xx[ 3 + nn*4+mm*4*Len ] est l'élément (nn,mm) de la matrice V-1 (poles en 1/eps qui vient du calcul après renormalisation -> divergence IR à compenser avec les PDF)

      res+=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      res-=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));

      // symetry case u ubar and ubar u
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);

      tu=1;
      res+=std::imag(TraceHSDiff(n,Xbet,chan,M2,xx2,tu)*ColinearDiff(n,chan,M2)*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      res-=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx2,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      
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
      
      ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
      
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);

      tu=1;
      res+=std::imag((TraceHSExpDiff(n,Xbet,chan,M2,xx2,tu)+(1.+ColinearExpDiff(n,chan,M2))*TraceBornDiff(n,Xbet,chan,M2))*MellinPDFDiff(n,chan,M2)*GlobalDiff(x,sc,M2));
    }
    else if(mel==4)
    {
      std::complex<double> tmp=0.;
      double a=-1.;
      tmp=TraceBornDiff(n-1.,Xbet,chan,M2)*2.; // *2 because of qqb <-> qbq symmetry 
      std::complex<double> g_12=std::pow(n+a-1.,-2.*2);// fake pdf
      res=M_PI*imag(tmp*g_12*xPDFDiff(x,chan,sc,M2,F,2)*GlobalDiffPDF(x,a,sc,M2));
    }

  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  return res;
}
