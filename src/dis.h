#ifndef DIS_H
#define DIS_H

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
using namespace Rcpp;
//using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// ===================================================

/**
* Class with variables and functions related to the discritization of the state space of the HMDP. Only store the variables needed (rest in R).
* \author Reza Pourmoayed & Lars Relund.
*/
class DIS
{
public:  // methods

    /** Constructor. Store the parameters.
     * @param dlmParam A list created using \code{iniDLM} and \code{buildDLM} in R.
     */
   DIS(const List param) {
      
      List rParam(param); 
      int n, k;
      List sWtmp = as<List>(rParam["wIntervals"]); 
      tMax = sWtmp.size();
      for (int i=0;i<tMax;i++) {
         NumericMatrix tmp = as<NumericMatrix>(sWtmp[i]);
         n = tmp.nrow(), k = tmp.ncol();
         sW.push_back( arma::mat(tmp.begin(), n, k, true) );
      }
      
      List sSdtmp = as<List>(rParam["sDIntervals"]); 
      for (int i=0;i<tMax;i++) {
         NumericMatrix tmp = as<NumericMatrix>(sSdtmp[i]);
         n = tmp.nrow(); k = tmp.ncol();
         sSd.push_back( arma::mat(tmp.begin(), n, k, true) );
      }
      List sGtmp = as<List>(rParam["gIntervals"]); 
      rations = sGtmp.size();
      for (int i=0;i<rations;i++) {
         NumericMatrix tmp = as<NumericMatrix>(sGtmp[i]);
         n = tmp.nrow(); k = tmp.ncol();
         sG.push_back( arma::mat(tmp.begin(), n, k, true) );
      }                        
   }
         
public:
 vector<arma::mat> sW;   // A list containing the matrices of of the possible weight intervals at the stages of the HMDP  
 vector<arma::mat> sG;   // A list containing the matrices of of the possible growth rates intervals at the stages of the HMDP  
 vector<arma::mat> sSd;  // A list containing the matrices of of the possible standard deviation intervals at the stages of the HMDP  
   
private:   
   int tMax;
   int rations;
};


#endif