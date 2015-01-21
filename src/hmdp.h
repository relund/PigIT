
#ifndef HMDP_HPP
#define HMDP_HPP

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "binaryMDPWriter.h"
#include "dlm.h"
#include "dglm.h"
#include "dis.h"
#include "cumNorm.h"
#include "time.h"
#include "limits"

using namespace Rcpp;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do

// [[Rcpp::depends(RcppArmadillo)]]

// ===================================================

/**
* Class for building a 3 level HMDP modelling a pen with info about weight and feeding.
* 
* It store all the info needed for writing the HMDP to the binary files.
* 
* @author Lars Relund (lars@relund.dk)
*/
class HMDP
{
public:  // methods

    /** Constructor. Store the parameters.
     * 
     * @param filePrefix Prefix used by the binary files storing the MDP model.
     * @param param Model parameters a list created using \code{setParameters} in R.
     * @param dlms A list a long as the number of rations where each entry is a DLM created using \code{iniDLM} and 
     *   \code{buildDLM} in R.
     * @param dglmParam Model parameters a list created using \code{iniDGLM} in R.
     */
   HMDP(const string filePrefix, const List param, const List dlms, const List dglmParam);
   
   
   /** Build the HMDP (to binary files). Use "shared linking". 
    *  
    *  Build a 3 level HMDP saved in binary files by using the 
    *  binaryMDPWriter (c++ version). The MDP specified at 3. level is generated for
    *  each ration for a dummy state at the first stage at level 2. .
    *  
    *  @return Build log (string)
    *  @author Lars Relund \email{lars@@relund.dk}
    *  @export
    */
   SEXP BuildHMDP();
   
   
   /** Count the number of states in the HMDP */
   int countStatesHMDP() {
      int ph, y, iRa;
      y=0;
      for(ph=0;ph<phases;ph++){
         if (ph==0) y=1; else y=y+IdCountL1(ph-1);
         for (iRa=0;iRa<rations;iRa++){  
               if( (ph==0) & (iRa!=iniFeedId) ) continue;
               y=y+IdCountL2(ph);
         }
      }
      return y + IdCountL1(phases-1) +1;
   }
   
   
private:

   /** Compare two vectors and return an error if not equal. */
   template <typename T>
   bool inline compareVec(const vector<T>& v1, const vector<T>& v2) {
      if (v1.size()!=v2.size()) {
         Rcout << "Vectors do not have the same size!" << endl;
         Rcout << "v1: " << vec2String<T>(v1) << " v2: " << vec2String<T>(v2) << endl;
         return false;
      }
      for (int i=0;i<v1.size(); i++) {
         if(!Equal(v1[i],v2[i])) {
            return false;
         }
      }
      return true;
   }
   
   /** Compare all input vectors with tmp vectors*/
   void compareAllVec(string state, string action, const vector<int> & idxTmp, const vector<int> & scpTmp, const vector<flt> prTmp) {
      bool error = false;
      if (!compareVec<flt>(prTmp,pr)) error = true;
      if (!compareVec<int>(idxTmp,index)) error = true;
      if (!compareVec<int>(scpTmp,scope)) error = true;
      if (error) {
         Rcout << "Error when comparing vectors in state " << state << " and action " << action << ":\n";
         Rcout << "idxTmp: " << vec2String<int>(idxTmp) << endl << " index: " << vec2String<int>(index) << endl;
         Rcout << "scpTmp: " << vec2String<int>(scpTmp) << endl << " scope: " << vec2String<int>(scope) << endl;
         Rcout << "prTmp: " << vec2String<flt>(prTmp) << endl << " pr: " << vec2String<flt>(pr) << endl;
         Rcpp::stop("");
      }
   }

   /** Allocate memory for matrices used given threshold actions.
    *  
    *  @author Lars Relund \email{lars@@relund.dk}
    */
   void AllocateMemTH();
   

   /** Allocate memory for matrices used given cull actions. */
   void AllocateMemCull();


   /** Calculate and fill arrays with rewards and trans pr. */
   void Preprocess();


   /** Create the process at level 1. */ 
   void BuildL1Process();
   
   
   /** Create the stages at level 1. */
   void BuildL1Phase();
   
   
   /** Create last stage of level 1 (dummy stage). */
   void BuildL1StageLast();
   
   
   /** Create the process at level 2 for a specific ration and phase when using threshold actions. 
    * 
    * @param iRation Ration index under consideration.
    */
   void BuildL2ProcessTH(int phase,int iRation);
   

   /** Create the process at level 2 for a specific ration and phase when using cull actions. 
    * 
    * @param iRation Ration index under consideration.
    */
   void BuildL2ProcessCull(int phase,int iRation);
    
    
   /** Calculate the initial transition probabilities and weights for the first week. 
    * 
    * Set the class vectors scope, pr and index based on global vector iniDist.
    */  
   void WeightTransPrIni();
   
   
   /** Calculate the weight values under threshold actions and Pr( n_t+1| n_t, t, iTH ) based on simulation .
   * 
   *  Values are stored in the 7-dim vector \var(weightTH[t][iTH][iR][iSW][iSG][iSSd][n]) and 6-dim vector \var(PrN[t][iTH][iSWt][iSSdt][nt][n]).
   */
   void CalcWeightsPrNTH();
   
   
   /** Calculate the probabilities Pr( n_t+1| n_t, t, iTH ).
    * 
    *  Values are stored in the 6-dim vector \var(PrN[t][iTH][iSWt][iSSdt][nt][n]).
    */
   void CalcPrN();
   
   
   /** Calculate the probabilities Pr(weight <= threshold) (unselected pen) used when we calculate Pr(n_t+1 | n_t).
    * 
    * Values are stored in the 3-dim vector \var(successPr[t][iTH][iSWt][iSSd]).
    */
   void CalcSuccessPr();
   
   
  /** Calculate the probabilities Pr( n_t+1| n_t, t, iTH ) based on simulation given a threshold action at level 2. 
    * 
    * @param iTH Index of action (entry in threshold vector).
    * @param nt Number of pigs in the pen at time t.
    */
   void CalcPrNSim(int & t, int & iTH, int & nt, int & iSWt, int & iSSdt, arma::mat & oWM, arma::mat & lWM, arma::mat & cWM);     
   

   /** Set the weights and calculate the transition probabilities of a threshold action at level 2. 
    * 
    * Set the class vectors weights, scope, pr and index. 
    * 
    * @param iRation Ration index under consideration. 
    * @param t Week number.
    * @param iTH Index of action (entry in threshold vector).
    * @param nt Number of pigs in the pen at time t.
    * @param iSWt Index of state variable for weight at time t.
    * @param iSGt Index of state variable for growth at time t.
    * @param iSSdt Index of state variable for std. dev at time t.
    */
   void WeightsTransPrTH(int & iRation, int & RSt, int & phase, int & t, int & iTH, int & nt, int & iSWt, int & iSGt, int & iSSdt);


   /** Set the weights and calculate the transition probabilities of a cull action at level 2. 
    * 
    * Set the class vectors weights, scope, pr and index. 
    * 
    * @param iRation Ration index under consideration. 
    * @param t Week number.
    * @param q Number of pigs to cull.
    * @param nt Number of pigs in the pen at time t.
    * @param iSWt Index of state variable for weight at time t.
    * @param iSGt Index of state variable for growth at time t.
    * @param iSSdt Index of state variable for std. dev at time t.
    */
   void WeightsTransPrCull(int & iRation, int & RSt, int & phase, int & t, int & iTH, int & nt, int & iSWt, int & iSGt, int & iSSdt);


  /** Sample the pigs in the unselected pen given a mean and variance given class members sDMeasure and samples. The samples are stored in 3 matrices.
    * 
    * @param mean Mean weight in the pen.
    * @param sd Std. dev. in the pen.
    * @param oWM arma matrix to store the sorted observed live weights for all the samples.
    * @param lWM arma matrix to store the true live weights for all the samples.
    * @param cWM arma matrix to store the carcass weights for all the samples.
    * 
    * @author Lars Relund \email{lars@@relund.dk}
    */
   void SimulatePigs(const double & mean, const double & sd, arma::mat & oWM, arma::mat & lWM, arma::mat & cWM);
   
   
  /** Calculate the reward of a threshold action at level 2. 
    * 
    * @param iTH Index of action (entry in threshold vector).
    * @param t Week number (stage).
    * @param iRt Ration index under consideration at time t. 
    * @param tRt Time used ration at time t.
    * @param nt Number of pigs in the pen at time t.
    * @param iSWt Index of state variable for weight at time t.
    * @param iSGt Index of state variable for growth at time t.
    * @param iSSdt Index of state variable for std. dev at time t.
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
    */
   double RewardTH(int & iRt, int & t, int & iTH, int & nt, int & iSWt, int & iSGt, int & iSSdt,
                        arma::mat & oWM, arma::mat & lWM, arma::mat & cWM); 
                        

  /** Calculate the weights of a cull action at level 2. 
    * 
    * @param cull number of cull pigs.
    * @param t Week number (stage).
    * @param iRt Ration index under consideration at time t. 
    * @param tRt Time used ration at time t.
    * @param nt Number of pigs in the pen at time t.
    * @param iSWt Index of state variable for weight at time t.
    * @param iSGt Index of state variable for growth at time t.
    * @param iSSdt Index of state variable for std. dev at time t.
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
    */
   double RewardCull(int & iRt, int & t, int & cull, int & nt, int & iSWt, int & iSGt, int & iSSdt,
                        arma::mat & oWM, arma::mat & lWM, arma::mat & cWM);
                        
   

   /** Calculate the transition probabilities of a terminate action at level 2. 
    * 
    * @param phase The phase number in the current stage
    * Set the class vectors scope, pr and index.
    */
   void WeightsTransPrTermTH(int & phase, int & iRation, int & t, int & n, int & iSW, int & iSG, int & iSSd) ;

   
   /** Calculate the transition probabilities of a change ration action at level 2 and set the weights to zero. 
    * 
    * Set the class vectors scope, pr and index.
    * 
    * @param iRation Ration index under consideration. 
    * @param t Week number.
    * @param iA Index of action (entry in threshold vector).
    * @param n Number of pigs in the pen at time t.
    * @param iSW Index of state variable for weight at time t.
    * @param iSSd Index of state variable for std. dev at time t.
    */
   void WeightsTransPrChange(int & iRationt, int & t, int & nt, int & iSWt, int & iSSdt);
   
      
   
   /** Calculate the log transformed probabilities Pr(m_t+1 | m_t, t).
    * 
    *  Values are stored in the 7-dim vector \var(prM[t][iRS][iSWt][iSGt][iSW][iSG][iR]).
    */
   void CalcTransPrM();
   
   
     /** Calculate the log transformed probabilities Pr(var_t+1 | var_t, t, n).
    * 
    *  Values are stored in the 3-dim vector \var(prSd[t][iSSdt][iSSd]).
    */
   void CalcTransPrSd();
      
      
     /** Calculate the  weight values under action Marketing.
    * 
    *  Values are stored in the 7-dim vector \var(weightTH[t][iTH][iR][iSW][iSG][iSSd][n] ).
    */
   void CalcWeightsCull();
        
   
   /** Build the map \var{mapL1}. That is, the map to identify state id at level 1 given string (t,n,iSW,iSSd).
    */ 
   void BuildMap(int phase);
   
    
    /** Build the map \var{mapL2}. That is, the map to identify state id at level 2 given current week to spped up the codes.
    */ 
   void BuildMapL2Vector(int week, int phase);
   
   
    /** Count the number of states in level 1. 
    */
    int IdCountL1(int phase);
    
     
    /** Count the number of states in level 2. 
    */
    int IdCountL2(int phase);
    
   
   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b) {
      std::ostringstream s;
      s << "(" << a << "," << b << ")";
      return s.str();
   }
   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << ")";
      return s.str();
   }
   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c, const int & d) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << "," << d << ")";
      return s.str();
   }
   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << "," << d << "," << e << ")";
      return s.str();
   }
   
   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e, const int & f) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << ")";
      return s.str();
   }

   /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e, const int & f, const int & g) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << "," << g << ")";
      return s.str();
   }
   
         /** Convert integers into a string. */
   string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e, const int & f, const int & g, const int & h) {
      std::ostringstream s;
      s << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << "," << g << "," << h << ")";
      return s.str();
   }


/** Compute the carcass price of one pigs based on given weight (Anders version) . 
    * @param  sWeight True carcass weight.
    * 
    * @return price of one pig with true carcass weight sWeight
    * 
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
    */    
   double PriceCarcass(double sWeight);

  
  /** Compute the lean meat price of one pigs based on given weight (Anders version) . 
    * @param sWeight True carcass weight.
    * @param tWeight True weight.
    * @param t Current week.
    * 
    * @return leannes price of one pig.
    * 
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
    */    
    double PriceLeanness(double sWeight, double tWeight, int t);
  
  
  
  /** Compute the feed cost of one pigs based on given weight . 
    * @param  tWeight True  weight.
    * @param  meanGrowth True  growth.
    * @param  iRt Index of feed-mix.
    * 
    * @return feed cost of one pig with true weight tWeight
    * 
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
    */
    double FeedCost(double tWeight, double meanGrowth, int iRt, int t);
    
    
//----------------------------------------------------------------------------------------------------------------------------------


      
private:   // variables

   static const double ZERO;  // trans pr below are considered as ZERO

   int phases;
   int pigs;
   int rations;
   int iniFeedId;
   int tMax, tStartMarketing;
   double CovWG;                    
   double convRate;               
   double convRateSd;             
   double avgGRate;                 
   double avgLeanP;               
   double avgInsWeight;           
   arma::vec rationCost;          
   arma::vec iniDist;           
   arma::vec priorGrowth;
   arma::vec rGIdx;             
   arma::vec tH;        // thresholds (kg)
   arma::vec minPhaseT; // minimum phases length 
   bool cullActions; // true if use cull actions instead of threshold actions
   int samples;      // number of samples when simulate   
   bool check;       // Do various checks e.g. if trans pr sum to one.
   double sDMeasure; // sd on measurements (e.g. video)
   int removeA;      // # of actions removed
   flt pigletCost;   // cost of a piglet
   
   int sizeTH;
   int sizeSW;
   int sizeSG;
   int sizeSSd;
   vector<DLM> dlm;     // DLMs used (one for each ration)
   DGLM dglm;           // DGLM model used 
   DIS dis;             // An object for the discritization class (DIS) 
   
   // variables used when build process
   vector<int> scope;
   vector<int> index;
   vector<flt> pr;
   vector<flt> weights;
   
   map<string,int> mapL1;   // map to identify state id at level 1 given string (t,n,iSW,iSSd) 
   map<string,int> mapR;    // find unique id for a state at level 2 (whole process)
   vector<vector<vector< vector< vector<int> > > > > mapL2Vector; // mapL2Vector[RS][n][iSW][iSG][iSSd] Vector of mapL2 instead of strings to speed up the codes in the second level 
   
   vector< vector< vector < vector < vector < vector < vector<double> > > > > > > prM; // prM[t][iRS][iSWt][iSGt][iSW][iSG][iR]  log probability (m_t+1|m_t, t).
   vector< vector< vector<double> > > prSd; // prSd[t][iSSdt][iSSd] log probability (var_t+1|var_t, t, n)
   vector< vector< vector< vector < vector <vector <vector <double> > > > > > > weightCull; // weightTH[t][iTH][iR][iSW][iSG][iSSd][n]  weight values for marketing decisions 
   
   vector<vector< vector< vector<flt> > > > successPr; // successPr[t][iTH][iSWt][iSSd]: probablity of success (used when calc trans pr for (n_t+1|n_t))   
   vector< vector<vector<vector< vector< vector<flt> > > > > > PrN; // PrN[t][iTH][iSWt][iSSdt][nt][n] log probabilities Pr( n_t+1| n_t, t, iTH ).
   vector< vector< vector< vector < vector <vector <vector <flt> > > > > > > weightTH; // weightTH[t][iTH][iR][iSW][iSG][iSSd][n]  weight values for threshold actions 
   
   string label;
   binaryMDPWriter w;
   ostringstream s;         // stream to write labels
   
   TimeMan cpuTime;
};


#endif