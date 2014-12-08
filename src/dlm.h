#ifndef DLM_H
#define DLM_H

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "cumNorm.h"
using namespace Rcpp;
//using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// ===================================================

/**
* Class with variables and functions related to a DLM. Only store the variables needed (rest in R).
* \author Lars Relund.
*/
class DLM
{
public:  // methods

    /** Constructor. Store the parameters.
     * @param dlmParam A list created using \code{iniDLM} and \code{buildDLM} in R.
     */
   DLM(const List dlmParam) {
      List rParam(dlmParam); 
      
      NumericMatrix tmp = as<NumericMatrix>(rParam["GG"]); 
      int n = tmp.nrow(), k = tmp.ncol();
      arma::mat tmp1(tmp.begin(), n, k, false);       // copy memory
      GG = tmp1;
      
      k1 = as<arma::vec>(dlmParam["X"]);
      k2 = as<double>(dlmParam["k2"]); 
      
      List cov = as<List>(rParam["L"]); 
      tMax = cov.size();
      for (int i=0;i<tMax;i++) {
         tmp = as<NumericMatrix>(cov[i]);
         n = tmp.nrow(); k = tmp.ncol();
         CovMat.push_back( arma::mat(tmp.begin(), n, k, true) );
      }
   }
   
   /** Transition probability, i.e. bivariate cumulative normal for two level. */  //Two levels
   double transPr(int week, int rationLength, arma::vec lower, arma::vec upper, arma::vec mt) { 
      //DBG3(endl << "    GG:" << GG << " mt:" << mt << endl)
      arma::vec mean= GG * mt;
      arma::mat sigma(2,2);
      sigma(0,0)=CovMat[week-1](0,0); sigma(0,1)=sigma(1,0)=CovMat[rationLength](0,1); sigma(1,1)=CovMat[rationLength](1,1);  
      DBG3(endl << "    pNorm2D: l=" << lower << " u=" << upper << "m=" << mean << " cov=" << sigma << " pr=" << pNorm2D_arma(lower, upper, mean, sigma) << " ")
      return pNorm2D_arma(lower, upper, mean, sigma);             
   }
   
      /** Transition probability, i.e. bivariate cumulative normal for three level. */ //Three levels 
   double transPr1(int week, int rationStart, arma::vec lower, arma::vec upper, arma::vec mt) { 
      //DBG3(endl << "    GG:" << GG << " mt:" << mt << endl)
      arma::vec mean = GG * mt;
      
      //DBG4("start: "<<endl)
      //DBG4("mt: "<<mt<<endl)
      //DBG4("mean: "<<mean<<endl)
      //DBG4("end: "<<endl)
      
      arma::mat sigma(2,2);       
       //mean[0] = mt[0] + mt[1];
       //mean[1] = mt[1];
       sigma(0,0)=CovMat[week-1](0,0); sigma(0,1)=sigma(1,0)=CovMat[week-rationStart](0,1); sigma(1,1)=CovMat[week-rationStart](1,1);  
       DBG3(endl << "    pNorm2D: l=" << lower << " u=" << upper << "m=" << mean << " cov=" << sigma << " pr=" << pNorm2D_arma(lower, upper, mean, sigma) << " ")
       return pNorm2D_arma(lower, upper, mean, sigma);  
   }
   
      /** Log transformed transition probability, i.e. bivariate cumulative normal. */
   double logTransPr(int week, int rationLength, arma::vec lower, arma::vec upper, arma::vec mt) { 
      return log( transPr(week, rationLength, lower, upper, mt) );
   }
   
   
        /** Log transformed transition probability, i.e. bivariate cumulative normal. */
   double logTransPr1(int week, int rationStart, arma::vec lower, arma::vec upper, arma::vec mt) { 
      return log( transPr1(week, rationStart, lower, upper, mt) );
   }
   
   
public:

   arma::vec k1;
   double k2;
   
private:   
   int tMax;
   arma::mat GG;
   vector<arma::mat> CovMat;   // covariance matrices of the transition probabilities. Note CovMat for week t+1 is stored in CovMat[week]!  

};


#endif