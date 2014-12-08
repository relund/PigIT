#ifndef DGLM_H
#define DGLM_H

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "cumNorm.h"
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// ===================================================

/**
* Class with variables and functions related to a DGLM. Only store the variables needed (rest in R).
* \author Lars Relund.
*/
class DGLM
{
public:  // methods

    /** Constructor. Store the parameters.
     * @param dglmParam A list created using \code{iniDGLM}.
     */
   DGLM(const List dglmParam) {
      List rParam(dglmParam); 
      beta0 = as<double>(rParam["beta0"]);
      numSample = as<int>(rParam["nf"]);
   }
   
   /** Transition probability, i.e. inverse gamma cdf. 
    *
    * @param t Week (time t).
    * @param n Number of pigs at time t.
    * @param lower Lower bound on variance at time t+1.
    * @param upper Upper bound on variance at time t+1.
    * @param var Estimate of variance at time t.
    */
   double logTransPr(int t, double lower, double upper, double var) { //SOLVED[Reza] : Based on the formulatiom for this probability in the paper, I changed "n" to "nf" (nf is the sample size). 

      double G;
      double probSd, xUpper, xLower; 
      double a, s, alpha, gamma, beta;
      if( (t==1) || (t==2) ){
         G = 1.8;
      }else{
         G = pow( (double) (t+1)/(t),2);
      }       
      double oShape = (double) (numSample-1)/(2);   //shape parameter of observation distribution
      double iShape = (double) (numSample-3)/(numSample-5); //("shape parameter of prior at t=1"): c_1 in the paper     
      
      a = (double) (G * var * ( iShape + oShape*t ) ) / ( iShape + oShape*(t+1) ) ; // location
      s = a; // scale
      alpha = oShape; // shape 1
      gamma = iShape + oShape*t +1; // shape 2
      beta =1 ;  // Weibul parameter
      xUpper = (double) (1)/( 1 + pow ( (double)(upper-a)/(s),-beta ) );
      xLower = (double) (1)/( 1 + pow ( (double)(lower-a)/(s),-beta ) );
      
      probSd=log( R::pbeta(xUpper, alpha, gamma,1, 0) - R::pbeta(xLower, alpha, gamma,1, 0) );
      if( ( R::pbeta(xUpper, alpha, gamma,1, 0) - R::pbeta(xLower, alpha, gamma,1, 0) )<0 ) DBG4("error_minus"<<endl)
      //if(probSd!=probSd)  DBG4(endl << " Error" << " lower=" << lower << " upper=" << upper << " centerp=" << var <<" t: "<<t<< endl)
      return (probSd);  // Rf_pgamma(q, shape, scale, lower.tail, log.p)
      }
       
public: 
   int numSample;


private:   
   double beta0;
};


#endif