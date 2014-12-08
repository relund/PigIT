#include "hmdp.h"

// ===================================================

const double HMDP::ZERO = 4.656613e-10;   // equal 1/2147483647 (limit of int since trans pr is stored as int when build the hgf) 

// ===================================================

HMDP::HMDP(const string prefix, const List param, const List dlms, const List dglmParam) : dglm(dglmParam), dis(param), w(prefix) {
   List rParam(param);       // Get parameters in params
   phases = as<int>(rParam["phases"]);
   pigs = as<int>(rParam["pigs"]);
   tMax = as<int>(rParam["tMax"]);
   tStartMarketing = as<int>(rParam["tStartMarketing"]);               
   rations = as<int>(rParam["rations"]);   
   iniFeedId = as<int>(rParam["iniFeedMix"])-1;
   minPhaseT = as<arma::vec>(rParam["minPhaseT"]); 
   CovWG = as<double>(rParam["CovWG"]);                  
   convRate = as<double>(rParam["convRate"]);            
   convRateSd = as<double>(rParam["convRateSd"]);        
   avgGRate = as<double>(rParam["avgGRate"]);            
   avgLeanP = as<double>(rParam["avgLeanP"]);            
   avgInsWeight = as<double>(rParam["avgInsWeight"]);    
   rationCost = as<arma::vec>(rParam["rationCost"]);      
   iniDist = as<arma::vec>(rParam["iniDist"]);         
   priorGrowth = as<arma::vec>(rParam["priorGrowth"]);
   rGIdx = as<arma::vec>(rParam["rGIdx"]);      
   tH = as<arma::vec>(rParam["thresholds"]); 
   cullActions = as<bool>(rParam["cullActions"]);
   samples = as<int>(rParam["samples"]);
   check = as<bool>(rParam["check"]);   
   sDMeasure = as<double>(rParam["sDMeasure"]); 
   removeA = 0;
   pigletCost = 375;
   
   sizeTH =tH.size();  
   sizeSW =dis.sW[0].n_rows; 
   sizeSG =dis.sG[0].n_rows;
   sizeSSd =dis.sSd[0].n_rows;
   
   List rDlms(dlms); 
   for (int i=0; i<rations; i++) {
      dlm.push_back( DLM( as<List>(rDlms[i])) );
   }
  
   // matrices for filling the rewards and transition probabilities before running the HMDP:
   prM = vector< vector< vector < vector < vector < vector < vector<double> > > > > > > (tMax, 
      vector< vector < vector < vector < vector < vector<double> > > > > >(tMax, 
      vector< vector < vector < vector < vector<double> > > > >(sizeSW, 
      vector< vector < vector < vector<double> > > >(sizeSG, 
      vector< vector < vector<double> > >(sizeSW,
      vector< vector<double> >(sizeSG,
      vector <double>(rations) ) ) ) ) ) );  // prM[t][iRS][iSWt][iSGt][iSW][iSG][iR]
   prSd = vector< vector< vector<double> > > (tMax, 
      vector< vector<double> >(sizeSSd, 
      vector<double>(sizeSSd) ) ); //prSd[t][iSSdt][iSSd]
   mapL2Vector = vector<vector< vector< vector< vector<int> > > > >(tMax, 
      vector<vector< vector< vector<int> > > >(pigs+1, 
      vector<vector< vector<int> > >(sizeSW,
      vector<vector<int> > (sizeSG, vector <int>(sizeSSd) ) ) ) ); //mapL2Vector[RS][n][iSW][iSG][iSSd]
}

// ===================================================

void HMDP::AllocateMemTH() {
   // matrices for filling the rewards and transition probabilities before running the HMDP:
   successPr = vector< vector< vector< vector<double> > > > (tMax, 
      vector< vector< vector<double> > >(sizeTH, 
      vector< vector<double> >(sizeSW,
      vector<double>(sizeSSd) ) ) ); //successPr[t][iTH][iSWt][iSSd]
   PrN= vector< vector<vector<vector< vector< vector<double> > > > > > (tMax, 
      vector< vector< vector<vector< vector<double> > > > > (sizeTH,
      vector< vector< vector<vector<double> > > >(sizeSW, vector< vector<vector<double> > >(sizeSSd,
      vector< vector<double> >(pigs+1, vector<double>(pigs+1) ) ) ) ) ); //PrN[t][iTH][iSWt][iSSdt][nt][n]   
   weightTH = vector< vector< vector< vector < vector <vector <vector <double> > > > > > > (tMax+1, 
      vector< vector< vector < vector <vector <vector <double> > > > > >(sizeTH,
      vector< vector < vector <vector <vector <double> > > > >(rations,
      vector < vector <vector <vector <double> > > >(sizeSW, 
      vector <vector <vector <double> > >(sizeSG,
      vector <vector <double> >(sizeSSd,
      vector <double>(pigs+1) ) )  ) ) )  ); //weightTH[t][iTH][iR][iSW][iSG][iSSd][n] 
}

// ===================================================

void HMDP::AllocateMemCull() {
   // matrices for filling the rewards and transition probabilities before running the HMDP:   
   weightCull = vector< vector< vector< vector < vector <vector <vector <double> > > > > > > (tMax+1, 
      vector< vector< vector < vector <vector <vector <double> > > > > >(pigs+1,
      vector< vector < vector <vector <vector <double> > > > >(rations,
      vector < vector <vector <vector <double> > > >(sizeSW, 
      vector <vector <vector <double> > >(sizeSG,
      vector <vector <double> >(sizeSSd,
      vector <double>(pigs+1) ) )  ) ) )  ); //weightTH[t][cull][iR][iSW][iSG][iSSd][n] 
}

// ===================================================

void HMDP::Preprocess() {
   if (cullActions) {
      Rcout << "Build the HMDP using cull actions ... \n\nStart preprocessing ...\n";
      AllocateMemCull();
      CalcWeightsCull();  
   }
   else {
      Rcout << "Build the HMDP using threshold actions ... \n\nStart preprocessing ...\n";
      AllocateMemTH();
      //CalcSuccessPr(); CalcPrN();  // Only if don't use sim
      CalcWeightsPrNTH();
   }
   CalcTransPrM();   
   CalcTransPrSd();
   Rcout << "... finished preprocessing.\n";
}

// ===================================================

SEXP HMDP::BuildHMDP() {
   Preprocess();
   
   Rcout << "Start writing to binary files ... \n";
   w.SetWeight("Time");
   w.SetWeight("Reward");
   
   w.Process();    // level 0 (founder)
   w.Stage();
   w.State("Dummy");
   WeightTransPrIni();  //calculate the initial transition probabilities
   w.Action(scope, index, pr, weights, "Dummy", false); 
      BuildL1Process();
   w.EndAction();
   w.EndState();
   w.EndStage();
   w.EndProcess();  // end level 1 (founder)
   w.CloseWriter();
   Rcout << "... finished writing to binary files.\n";
   Rcout << "Removed " << removeA << " actions.\n";
   return wrap(w.log.str());
}

// ===================================================


void HMDP::BuildL1Process() {
   w.Process(); // level 1
      BuildL1Phase();
      BuildL1StageLast();
   w.EndProcess(); // end level
}


// ===================================================


void HMDP::BuildL1Phase() {
   int n, iSW, iSSd, iA, iR, RS, phase;
      for(phase=0;phase<phases;phase++){         
         if ( phase<(phases-1) ) BuildMap(phase+1);
         w.Stage();
            // Dummy state to define the subprocesses
            w.State("Dummy");
               for (iR=0;iR<rations;iR++) { // rations
                  if( (phase==0) & (iR!=iniFeedId) ) continue;
                  scope.assign(1,2); index.assign(1,0); pr.assign(1,1); 
                  weights.assign(2,0); 
                  label = "iR" + ToString<int>(iR); 
                  w.Action(scope, index, pr, weights, label, false);
                     if (cullActions) BuildL2ProcessCull(phase,iR);
                     else BuildL2ProcessTH(phase,iR);
                  w.EndAction();
               }      
            w.EndState();
              
            for (iR=0;iR<rations;iR++) {  // rations
              if( (phase==0) & (iR!=iniFeedId) ) continue;
              for (RS=minPhaseT[phase];RS<tMax;RS++) { // possible time for changing the rations
              if( (phase==0) & (RS!=1) ) continue;
                  for (n=1;n<=pigs;n++) { // pigs in the pen
                     if( (RS<=tStartMarketing) & (n!=pigs) ) continue;
                     for (iSW=0; iSW<sizeSW; iSW++) {   // row number in sW (C style)
                        for (iSSd=0; iSSd<sizeSSd; iSSd++) {   // row number in sSd
                           label = getLabel(iR,RS,n,iSW,iSSd);
                           w.State(label);
                           for (iA = 0; iA<rations; iA++){
                              if ( (phase==0) & (iA!=iR) )  continue;   // only same ration for phase=0
                              if ( (phase!=0) & (iA==iR) )  continue;   // if same ration for phase>0 
                                 index.assign(1,mapR[getLabel(phase,iA,RS,n,iSW,rGIdx[iA],iSSd)]); scope.assign(1,3); pr.assign(1,1);  // sizeSG/2; instead of rGIdx[iA], index value for new ration is (sizeSG)/2 that always will be an integer  
                                 weights.assign(2,0);  
                                 label = "new fm = " + ToString<int>(iA); //
                                 w.Action(scope, index, pr, weights, label, true);      // action link to shared process
                           }  
                           w.EndState();
                        }
                     }
                  }
               }
            }    
            w.State("dummy");
               scope.assign(1,0); index.assign(1,0); pr.assign(1,1); weights.assign(2,0);
               w.Action(scope, index, pr, weights, "terminated", true);
            w.EndState();
         w.EndStage();
      } 
 }


// ===================================================


void HMDP::BuildL1StageLast() {
   w.Stage();
      w.State("dummy");
         scope.assign(1,0); index.assign(1,0); pr.assign(1,1); 
         weights.assign(2,0);  
         w.Action(scope, index, pr, weights, "terminated(s)", true);
      w.EndState();
   w.EndStage();
}

// ===================================================
void HMDP::BuildL2ProcessCull(int phase, int iRation) {
   int t, n, iSW, iSG, iSSd, RS, cull, idS;
   w.Process(); // level 2
      for(t = minPhaseT[phase]; t<=tMax; t++) {
         if(t != tMax) BuildMapL2Vector(t+1,phase);
         w.Stage();
             for(RS=minPhaseT[phase];RS<=t;RS++){
             if( (phase==0) & (RS!=1) ) continue;
             if(RS==tMax) continue;
               for (iSW=0; iSW<sizeSW; iSW++) {   
                  for (iSG=0; iSG<sizeSG; iSG++) {   
                     for (iSSd=0; iSSd<sizeSSd; iSSd++) {   
                       for (n=1;n<=pigs;n++) { // [Reza]: I changed 0 to 1. 
                         if ( (t<=tStartMarketing) & (n!=pigs) ) continue;     // currently no mortality -> always n=pigs in the pen before marketing
                         label = getLabel(n,iSW,iSG,iSSd,RS, phase, iRation, t);
                         idS = w.State(label);
                           if (t==RS) mapR[getLabel(phase,iRation,RS,n,iSW,iSG,iSSd)] = idS;
                            
                           if (t!=tMax) { 
                              cull = 0; // cull zero pigs action (cont.)
                              WeightsTransPrCull(iRation,RS,phase,t,cull,n,iSW,iSG,iSSd);
                              w.Action(scope, index, pr, weights, "cont.", true);
                              if (t>=tStartMarketing) {
                                 for (cull=1; cull<=n; cull++) {
                                    WeightsTransPrCull(iRation,RS,phase,t,cull,n,iSW,iSG,iSSd);
                                    label = "cull " + ToString<int>(cull);
                                    w.Action(scope, index, pr, weights, label, true);
                                 }
                              }
                           } 
                           else { // last stage
                              WeightsTransPrCull(iRation,RS,phase,t,n,n,iSW,iSG,iSSd);
                              w.Action(scope, index, pr, weights, "term.", true);
                           }
                           // change ration
                           if(  (3<=(t-RS)) & (t!=tMax) & (phase!=(phases-1)) ){ 
                               WeightsTransPrChange(iRation, t, n, iSW, iSSd);
                               w.Action(scope, index, pr, weights, "change fm", true);     
                           }   
                        w.EndState();
                     }
                  }
               }
            }
         }   
         w.EndStage();
      }
   w.EndProcess(); // end level        
}

// ===================================================
void HMDP::BuildL2ProcessTH(int phase, int iRation) { 
   int t, n, iSW, iSG, iSSd, RS, iTH, idS;
   //bool useA;  // use action
   w.Process(); // level 2
      for(t = minPhaseT[phase]; t<=tMax; t++) {
         if(t != tMax) BuildMapL2Vector(t+1,phase);
         w.Stage();
             for(RS=minPhaseT[phase];RS<=t;RS++){
               if( (phase==0) & (RS!=1) ) continue;
               if(RS==tMax) continue;
               for (iSW=0; iSW<sizeSW; iSW++) {   
                  for (iSG=0; iSG<sizeSG; iSG++) {   
                     for (iSSd=0; iSSd<sizeSSd; iSSd++) {   
                       for (n=1;n<=pigs;n++) { // [Reza]: I changed 0 to 1. 
                         if ( (t<=tStartMarketing) & (n!=pigs) ) continue;     // currently no mortality -> always n=pigs in the pen before marketing
                         label = getLabel(n,iSW,iSG,iSSd,RS, phase, iRation, t);
                         idS = w.State(label);
                           if (t==RS) mapR[getLabel(phase,iRation,RS,n,iSW,iSG,iSSd)] = idS;
                           
                           //We assume that iTH=0 is the term. action and that iTH=sizeTH-1 is the cont. action!! 
                           if (t!=tMax) { // if not last stage
                              // cont. action (iTH = sizeTH-1)
                              iTH = sizeTH-1;
                              WeightsTransPrTH(iRation,RS,phase,t,iTH,n,iSW,iSG,iSSd);
                              //CalcTransPrTH(iRation,RS,phase,t,iTH,n,iSW,iSG,iSSd);  //cout << getLabel(RS,t,iTH,n,iSW,iSG,iSSd) << "::" << flush; //cout << "scp:" << vec2String<int>(scope) << " idx:" << vec2String<int>(index) << " pr" << vec2String<flt>(pr) << " - " << flush;
                              //weights[0]=1; weights[1]=weightTH[t][iTH][iRation][iSW][iSG][iSSd][n];
                              label = ToString<double>(tH[iTH]);
                              w.Action(scope, index, pr, weights, label, true);
                              // true threshold actions
                              if (t>=tStartMarketing) {
                                 for (iTH=sizeTH-2; iTH>0; --iTH) { // scan thresholds backward     
                                    if (weightTH[t][iTH][iRation][iSW][iSG][iSSd][n]!=weightTH[t][iTH-1][iRation][iSW][iSG][iSSd][n]) { // weights are the same -> remove action
                                       WeightsTransPrTH(iRation,RS,phase,t,iTH,n,iSW,iSG,iSSd);
                                       //weights[0]=1; weights[1]=weightTH[t][iTH][iRation][iSW][iSG][iSSd][n];
                                       //if(  (weights[1])!=(weights[1])  ) DBG4("error_weight_mark"<<" t: "<<t<<" iSW: "<<iSW<<" iTH: "<<iTH<<endl)
                                       label = ToString<double>(tH[iTH]);
                                       w.Action(scope, index, pr, weights, label, true);
                                    } else removeA++;
                                 }
                                 // term. action (iTH = 0)
                                 iTH=0;
                                 WeightsTransPrTermTH(phase, iRation, t, n, iSW, iSG, iSSd);
                                 label = ToString<double>(tH[iTH]);
                                 w.Action(scope, index, pr, weights, label, true);
                              } 
                           } 
                           else { // last stage
                              WeightsTransPrTermTH(phase, iRation, t, n, iSW, iSG, iSSd);
                              w.Action(scope, index, pr, weights, "term.", true);
                           }
                           
                           // Decision at the end of phase for changing a ration
                           if(  (3<=(t-RS)) & (t!=tMax) & (phase!=(phases-1)) ){ 
                               WeightsTransPrChange(iRation, t, n, iSW, iSSd);
                               label= "change fm";
                               w.Action(scope, index, pr, weights, label, true);     
                           }          
                        w.EndState();
                     }
                  }
               }
            }
         }   
         w.EndStage();
      }
   w.EndProcess(); // end level        
}

// ===================================================

void HMDP::WeightsTransPrCull(int & iRation, int & RSt, int & phase, int & t, int & q, int & nt, int & iSWt, int & iSGt, int & iSSdt) {
   double pr4;
   int n, iSW, iSG, iSSd, RS, id;
   string str;
   RS=RSt;
      
   n=nt-q;  // n is the number of pigs after culling
   pr.clear(); index.clear(); scope.clear();
   if (q!=nt) {  // have not culled all pigs
      for (iSW=0; iSW<sizeSW; iSW++) {   
         for (iSG=0; iSG<sizeSG; iSG++) {   
            for (iSSd=0; iSSd<sizeSSd; iSSd++) {   
               id=mapL2Vector[RS][n][iSW][iSG][iSSd];
//               Rcout << id << " - ";
               pr4 = exp(prM[t][RS][iSWt][iSGt][iSW][iSG][iRation] + prSd[t][iSSdt][iSSd]);  
               //DBG4("info:" << getLabel(t,q,nt,n,iSW,iSG,iSSd,id) << " pr:" << ToString<flt>(prM[t][RS][iSWt][iSGt][iSW][iSG][iRation]) << "," << ToString<flt>(prSd[t][iSSdt][iSSd]) << "  ")
               if (pr4>ZERO) {
                  pr.push_back(pr4); index.push_back(id);
               }
            }
         }
      }
      weights[0] = 1; scope.assign(pr.size(),1); 
   } 
   else {  // cull all pigs
      pr4 = 1;
      if( phase!=(phases-1) ) id=mapL1["dummy"]; else id=0;
      pr.push_back(pr4); index.push_back(id); scope.push_back(0);
      weights[0] = 0; // cull instantly
   }
   weights[1] = weightCull[t][q][iRation][iSWt][iSGt][iSSdt][nt];  // set reward
   if (check) {
      arma::vec tmp(pr);
      if (!Equal(sum(tmp),1,1e-8)) {
         Rcout << "Warning sum pr!=1 in WeightsTransPrCull - diff = " << 1-sum(tmp) << " pigs = " << nt << " cull = " << q << " index:" << vec2String<int>(index) << " pr:" << vec2String<flt>(pr) << endl;
      }
      if (scope.size()==0) {
         str = "Error vectors empty in WeightsTransPrCul " + getLabel(t,q,nt,n,iSW,iSG,iSSd,id) + " pr: ("  + ToString<flt>(prM[t][RS][iSWt][iSGt][iSW][iSG][iRation]) + "," + ToString<flt>(prSd[t][iSSdt][iSSd]) + ")\n";
         Rcpp::stop(str);
      }
   }
}

// ===================================================

void HMDP::WeightsTransPrTH(int & iRation, int & RSt, int & phase, int & t, int & iTH, int & nt, int & iSWt, int & iSGt, int & iSSdt) {
   double pr1,pr4;
   int n, iSW, iSG, iSSd, RS, id;
   string str;
   RS=RSt;
      
   pr.clear(); index.clear(); scope.clear();
   for (n=1;n<=nt;n++) { 
      pr1 = PrN[t][iTH][iSWt][iSSdt][nt][n];
      if (exp(pr1)<ZERO) continue; 
      for (iSW=0; iSW<sizeSW; iSW++) {   
         for (iSG=0; iSG<sizeSG; iSG++) {   
            for (iSSd=0; iSSd<sizeSSd; iSSd++) {
               if (exp(prSd[t][iSSdt][iSSd])<ZERO) continue;
               id=mapL2Vector[RS][n][iSW][iSG][iSSd];
               pr4 = exp (pr1 + prM[t][RS][iSWt][iSGt][iSW][iSG][iRation] + prSd[t][iSSdt][iSSd]);  
               //DBG4("info:" << getLabel(t,iTH,nt,n,iSW,iSG,iSSd,id) << " pr:" << ToString<flt>(pr1) << "," << ToString<flt>(prM[t][RS][iSWt][iSGt][iSW][iSG][iRation]) << "," << ToString<flt>(prSd[t][iSSdt][iSSd]) << "  ")
               if (pr4>ZERO) {
                  pr.push_back(pr4);
                  index.push_back(id);
               }
            }
         }
      }
   }   
   scope.assign(pr.size(),1); 
   weights[0] = 1; // default time 
   weights[1] = weightTH[t][iTH][iRation][iSWt][iSGt][iSSdt][nt];  // set weight
   // zero pigs = trans to father
   n = 0;
   pr4=exp(PrN[t][iTH][iSWt][iSSdt][nt][n]);
   if (pr4>ZERO) {
      if( phase!=(phases-1) ) id=mapL1["dummy"];
      else id=0;
      pr.push_back(pr4);
      index.push_back(id);
      scope.push_back(0);
      weights[0] = 1-pr4; // set time if terminate with pr>0 
   } 
   if (check) {
      arma::vec tmp(pr);
      if (!Equal(sum(tmp),1,1e-7)) {
         Rcout << "Warning sum pr!=1 in WeightsTransPrTH - diff = " << 1-sum(tmp) << " pigs = " << nt << " iTh = " << iTH << " index:" << vec2String<int>(index) << " pr:" << vec2String<flt>(pr) << endl;
      }
      if (scope.size()==0) {
         str = "Error vectors empty in WeightsTransPrTH " + getLabel(t,iTH,nt,n,iSW,iSG,iSSd,id) + " pr: ("  + ToString<flt>(prM[t][RS][iSWt][iSGt][iSW][iSG][iRation]) + "," + ToString<flt>(prSd[t][iSSdt][iSSd]) + ")\n";
         Rcpp::stop(str);
      }
   }
}

// ===================================================

void HMDP::WeightsTransPrChange(int & iRationt, int & t, int & nt, int & iSWt, int & iSSdt) {
   int id, RS;
   RS=t;
   id=mapL1[getLabel(iRationt,RS,nt,iSWt,iSSdt)];
   scope.assign(1,0); 
   index.assign(1, id);
   pr.assign(1,1); 
   weights.assign(2,0); 
} 

// ===================================================

void HMDP::WeightTransPrIni() {   
   int iR, RS, n, iSW, iSSd;
   int id=1;
   int phaseI=0;
   double PrIni, PrMean, PrSd;
   pr.clear(); index.clear(); scope.clear();
    for (iR=0;iR<rations;iR++){
       if( (iR!=iniFeedId) ) continue;
       for (RS=minPhaseT[phaseI];RS<tMax;RS++){
          if( (RS!=1) ) continue;
           for (n=1;n<=pigs;n++){
              if( (RS<=tStartMarketing) & (n!=pigs) ) continue;
              for (iSW=0; iSW<sizeSW; iSW++){
                 for (iSSd=0; iSSd<sizeSSd; iSSd++){
                    PrMean=R::pnorm(dis.sW[0](iSW,2),iniDist[0],sqrt(iniDist[1]),1,0) - R::pnorm(dis.sW[0](iSW,1),iniDist[0],sqrt(iniDist[1]),1,0);
                    PrSd=exp( dglm.logTransPr( 1, pow(dis.sSd[0](iSSd,1),2), pow(dis.sSd[0](iSSd,2),2), iniDist[1] ) ); 
                    PrIni=PrMean*PrSd;
                    pr.push_back(PrIni);
                    index.push_back(id);
                    id++;
                 }
              }
           }
       }
    }
    scope.assign(pr.size(),2); 
    weights.assign(2,0); weights[1]=-pigletCost*pigs;
  } 

// ===================================================

void HMDP::BuildMap(int phase) {
   int n, iSW, iSSd, id, iR, RS;
   // level 1
   id=1;
   mapL1.clear();
      
   for (iR=0;iR<rations;iR++){    
       if( (phase==0) & (iR!=iniFeedId) ) continue;
       for (RS=minPhaseT[phase];RS<tMax;RS++) { 
         if( (phase==0) & (RS!=1) ) continue;         
         for (n=1;n<=pigs;n++){   
            if( (RS<=tStartMarketing) & (n!=pigs)) continue;
            for (iSW=0; iSW<sizeSW; iSW++){    
               for (iSSd=0; iSSd<sizeSSd; iSSd++){  
                  mapL1[getLabel(iR,RS,n,iSW,iSSd)] = id;
                  id++;
               }
            }
         }
      }
   }   
   mapL1["dummy"] = id;   
}


// =================================================== 

void HMDP::BuildMapL2Vector(int week, int phase) {   //  mapL2Vector[n][iSW][iSG][iSSd]
   int n, iSW, iSSd, iSG, RS, id;
   // level 2
   id=0;
   
      for(RS=minPhaseT[phase];RS<=week;RS++){
      if( (phase==0) & (RS!=1) ) continue;
      if(RS==tMax) continue;
      for (iSW=0; iSW<sizeSW; iSW++){
         for (iSG=0; iSG<sizeSG; iSG++){
            for (iSSd=0; iSSd<sizeSSd; iSSd++){
               for (n=1;n<=pigs;n++){
                  if ( (week<=tStartMarketing) & (n!=pigs) ) continue;
                   //mapL2Vector[RS][n][iSW][iSG][iSSd]=-1; 
                   mapL2Vector[RS][n][iSW][iSG][iSSd]= id;    
                   id++;
               }
            }
         }
      }
   }
}

// ===================================================

void HMDP::CalcWeightsCull(){   
   cpuTime.Reset(0); cpuTime.StartTime(0);
   int cull, iSW, iSSd, t, n, iR, iSG;
   double mean, sd;
   
   arma::mat oWM(samples,pigs), lWM(samples,pigs), cWM(samples,pigs);
   for(t = 1; t<tMax+1; ++t){
      for (iSW=0; iSW<sizeSW; iSW++){
         for (iSSd=0; iSSd<sizeSSd; iSSd++){
            mean = dis.sW[t-1](iSW,0);
            sd = dis.sSd[t-1](iSSd,0);
            SimulatePigs(mean, sd, oWM, lWM, cWM);   // store simulation in the matrices. Note use same sim for different th -> the rewards will tbe decreasing for increasing th
            for (iSG=0; iSG<sizeSG; iSG++){
               for (iR = 0; iR<rations; ++iR){
                  for(n=1;n<=pigs;n++) {
                     for (cull=0;cull<=n;cull++){
                        weightCull[t][cull][iR][iSW][iSG][iSSd][n] = RewardCull(iR,t,cull,n,iSW,iSG,iSSd, oWM, lWM, cWM);  // 2*t*iR*q*iSW*iSG*iSSd;
                     }
                  }
               }
            }
         }
      }
   }
   Rcout << "Time for calculating weightCull values (" << samples << " samples): " << cpuTime.StopAndGetTotalTimeDiff(0) << endl; 
}

// ===================================================

void HMDP::CalcTransPrM() {   // calc values prM[t][iRS][iSWt][iSGt][iSW][iSG][iR]
   cpuTime.Reset(0); cpuTime.StartTime(0);
   int t, iSWt, iSGt, iSW, iSG, iR, RS;
   arma::vec lower(2), upper(2), mt(2);
    for(t = 1; t<tMax; ++t) {
             for (iSWt=0; iSWt<sizeSW; iSWt++) {  
               for (iSGt=0; iSGt<sizeSG; iSGt++) {   
                 for(RS=1; RS<=t; RS++){
                    for (iR = 0; iR<rations; ++iR){
                      mt[0] = dis.sW[t-1](iSWt,0); mt[1] = dis.sG[iR](iSGt,0);
                      for (iSW=0; iSW<sizeSW; iSW++) {   
                        for (iSG=0; iSG<sizeSG; iSG++) {  
                         lower[0] = dis.sW[t](iSW,1); lower[1] = dis.sG[iR](iSG,1);  // lower at time t+1
                         upper[0] = dis.sW[t](iSW,2); upper[1] = dis.sG[iR](iSG,2);  // upper at time t+1
                         prM[t][RS][iSWt][iSGt][iSW][iSG][iR] = dlm[iR].logTransPr1(t, RS,lower, upper, mt);    // use log transf to avoid underflow
                     } 
                   }
                }
             }
          }
       }
    }
    Rcout << "Time for calculating prM: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl;
}

// ===================================================

void HMDP::CalcTransPrSd() {  // calc values for prSd[t][iSSdt][iSSd]
   cpuTime.Reset(0); cpuTime.StartTime(0);
   for(int t = 1; t<tMax; ++t){
      for (int iSSdt=0; iSSdt<sizeSSd; iSSdt++) {
         for (int iSSd=0; iSSd<sizeSSd; iSSd++) {
            prSd[t][iSSdt][iSSd] = dglm.logTransPr(t, pow(dis.sSd[t](iSSd,1),2), pow(dis.sSd[t](iSSd,2),2), pow(dis.sSd[t-1](iSSdt,0),2));
         }
      }   
   } 
   Rcout << "Time for calculating prSd: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl;
}  

// ===================================================

void HMDP::SimulatePigs(const double & mean, const double & sd, arma::mat & oWM, arma::mat & lWM, arma::mat & cWM) {
  oWM = arma::randn<arma::mat>(samples,pigs);   // use std. norm dist arma function
  oWM = oWM*sqrt(pow(sd,2) + pow(sDMeasure,2) ) + mean; // transform to general norm dist
  if (mean<50) oWM = abs(oWM); 
  oWM = arma::sort(oWM,"ascend",1);   // sort each row
  lWM = oWM - arma::randn<arma::mat>(samples,pigs)*sDMeasure; // NOTE TO DO: sdMeasure should actually be calculated from the V covariance matrix
  cWM = 0.84*lWM-5.89 + arma::randn<arma::mat>(samples,pigs)*convRateSd;  
}

// ===================================================

double HMDP::RewardCull(int & iRt, int & t, int & cull, int & nt, int & iSWt, int & iSGt, int & iSSdt,
                     arma::mat & oWM, arma::mat & lWM, arma::mat & cWM) {
   double simSum=0;
   double sumReveneue, sumCost;
   double meanGrowth = dis.sG[iRt](iSGt,0);
   int r,c;

   for (r=0; r<oWM.n_rows; r++) { // for each sample
      sumReveneue=0; sumCost=0; 
      for (c=0; c<nt-cull; c++) // for each pig 1,...,nt-q (keep the pigs)
         sumCost += FeedCost(lWM(r,c), meanGrowth, iRt, t);
      for(c=nt-cull; c<nt; c++)  // for each pig nt-q+1,...,nt (cull the q biggest)
         sumReveneue += PriceCarcass(cWM(r,c)) + PriceLeanness(cWM(r,c),lWM(r,c),t); 
      simSum += sumReveneue - sumCost;
   } 
   return(simSum/oWM.n_rows);
}  

// ===================================================

void HMDP::CalcWeightsPrNTH(){   
   cpuTime.Reset(0); cpuTime.StartTime(0);
   int iTH, iSW, iSSd, t, n, iR, iSG;
   double mean, sd;
   
   arma::mat oWM(samples,pigs), lWM(samples,pigs), cWM(samples,pigs);
   for(t = 1; t<tMax+1; ++t){
      for (iSW=0; iSW<sizeSW; iSW++){
         for (iSSd=0; iSSd<sizeSSd; iSSd++){
            mean = dis.sW[t-1](iSW,0);
            sd = dis.sSd[t-1](iSSd,0);
            SimulatePigs(mean, sd, oWM, lWM, cWM);   // store simulation in the matrices. Note use same sim for different th -> the rewards will tbe decreasing for increasing th
            for (iTH=0;iTH<sizeTH;iTH++){
               for(n=1;n<=pigs;n++) {
                  if (t<tMax) CalcPrNSim(t, iTH, n, iSW, iSSd, oWM, lWM, cWM); 
                  for (iR = 0; iR<rations; ++iR){
                     for (iSG=0; iSG<sizeSG; iSG++){
                        weightTH[t][iTH][iR][iSW][iSG][iSSd][n]= RewardTH(iR,t,iTH,n,iSW,iSG,iSSd, oWM, lWM, cWM);
                     }
                  }
               }
            }
         }
      }
   }
   Rcout << "Time for calculating weightTH values (" << samples << " samples): " << cpuTime.StopAndGetTotalTimeDiff(0) << endl; 
}

// ===================================================

void HMDP::CalcPrNSim(int & t, int & iTH, int & nt, int & iSWt, int & iSSdt, arma::mat & oWM, arma::mat & lWM, arma::mat & cWM) {
   double sum;
   double TH = tH[iTH];
   arma::vec kept(nt+1);  // kept[i] = # of samples where i pigs kept 
   
   kept.fill(0); 
   for (int r=0; r<oWM.n_rows; r++) { // for each sample
      sum = 0;
      for(int c=0; c<nt; c++) { // for each pig 1,...,nt
        if(oWM(r,c)<=TH) sum++;
      }
      kept[sum]++;
   }      
   kept = kept/oWM.n_rows;
   for(int n=0;n<=nt;n++){
      //Rcout << PrN[t][iTH][iSWt][iSSdt][nt][n] << ", ";
      PrN[t][iTH][iSWt][iSSdt][nt][n] = log(kept[n]);
      //if (!Equal(arma::sum(kept),1)) Rcout << getLabel(t,iTH,iSWt,iSSdt,nt,n) << ":" << kept << "   " ;
   }
}

// ===================================================

double HMDP::RewardTH(int & iRt, int & t, int & iTH, int & nt, int & iSWt, int & iSGt, int & iSSdt,
                     arma::mat & oWM, arma::mat & lWM, arma::mat & cWM) {
   double simSum=0;
   double sumReveneue, sumCost;
   double TH = tH[iTH];
   double meanGrowth = dis.sG[iRt](iSGt,0);
         
   for (int r=0; r<oWM.n_rows; r++) { // for each sample
      sumReveneue=0; sumCost=0; 
      for(int c=0; c<nt; c++) { // for each pig 1,...,nt
        if(oWM(r,c)>TH) {
           sumReveneue += PriceCarcass(cWM(r,c)) + PriceLeanness(cWM(r,c),lWM(r,c),t); // cull pig
        }
        else sumCost += FeedCost(lWM(r,c), meanGrowth, iRt, t);
      }
      simSum = simSum + sumReveneue - sumCost;
   }      
   return(simSum/oWM.n_rows);
}

// ===================================================

void HMDP::CalcPrN(){  
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iTH, iSWt, iSSdt, t, n, nt;
  string str;
  double sum;
  for(t = 1; t<tMax; ++t){
     for (iTH=0;iTH<sizeTH;iTH++){
        for (iSWt=0; iSWt<sizeSW; iSWt++){
           for (iSSdt=0; iSSdt<sizeSSd; iSSdt++){
              for(nt=1;nt<=pigs;nt++){
                 sum = 0;
                 for(n=0;n<=nt;n++){
                    PrN[t][iTH][iSWt][iSSdt][nt][n]=  // first calc for unselected pen
                       R::dbinom(n, (double)pigs, successPr[t][iTH][iSWt][iSSdt], 0/*no log transform*/);
                    sum+=PrN[t][iTH][iSWt][iSSdt][nt][n];
                 }
                 // next normalize
                 if (sum==0) {   // i.e. successPr[t][iTH][iSWt][iSSdt] = 1 and threshold so high that no pigs culled
                    PrN[t][iTH][iSWt][iSSdt][nt][nt] = 1;
                 } 
                 else for(n=0;n<=nt;n++) {
                    PrN[t][iTH][iSWt][iSSdt][nt][n] = PrN[t][iTH][iSWt][iSSdt][nt][n]/sum;
                 }
                 // log transform
                 for(n=0;n<=nt;n++) PrN[t][iTH][iSWt][iSSdt][nt][n] = log(PrN[t][iTH][iSWt][iSSdt][nt][n]);
              }
           }
        }
     }
  }
  Rcout << "Time for calculating PrN: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl; 
} 

// ===================================================

void HMDP::CalcSuccessPr() {  //successPr[t][iTH][iSWt][iSSd]
   cpuTime.Reset(0); cpuTime.StartTime(0);
   int iTH, iSWt, iSSdt, t;
   double sd; //successPr.clear();
   for(t = 1; t<tMax; ++t){
      for (iTH=0;iTH<sizeTH;iTH++) { 
         for (iSWt=0; iSWt<sizeSW; iSWt++) {   
            for (iSSdt=0; iSSdt<sizeSSd; iSSdt++) {   
               sd = sqrt(pow(dis.sSd[t-1](iSSdt,0),2) + pow(sDMeasure,2) ); 
               successPr[t][iTH][iSWt][iSSdt] = cumNorm1D( (tH[iTH]-dis.sW[t-1](iSWt,0))/sd );
               //R::pnorm(tH[iTH],dis.sW[t-1](iSWt,0),dis.sSd[t-1](iSSdt,0),1,0);
            }
         }
      }      
   }
   Rcout << "Time for calculating successPr: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl;   
}

// ===================================================

double HMDP::PriceCarcass(double sWeight){
   //Interval [0,60) :
   if(0<=sWeight && sWeight<60) return( sWeight*( 10.3 - 0.1*(70-60) - 0.2*(60-sWeight) ) );
   //Interval [60,70) :  
   if(60<=sWeight && sWeight<70) return( sWeight*( 10.3 - 0.1*(70-sWeight) ) );
   //Interval [70,86) : 
   if(70<=sWeight && sWeight<86) return( 10.3*sWeight);
   //Interval [86,95) :  
   if(86<=sWeight && sWeight<95) return( sWeight*( 10.3 - 0.1*(sWeight-86) ) );
   //Interval [95,100) :
   if(95<=sWeight && sWeight<100) return(9.3*sWeight);
   //Interval [100,Inf) : 
   if(100<=sWeight) return(9.1*sWeight);
   return(0);
}

// ===================================================

double HMDP::PriceLeanness(double sWeight, double tWeight, int t){
   double avgGrowth = (double) (tWeight - avgInsWeight)/(t);
   double lean = ( (double) (-7.5)/(7) )*(avgGrowth - avgGRate) + avgLeanP;
   return(sWeight*(lean - 61)*0.1 );
}

// ===================================================

double HMDP::FeedCost(double tWeight, double meanGrowth, int iRt, int t){
    
   //Cost based on daily feed intake and orginal parameters (version 1):
   /*double expGain = meanGrowth;
   double sum=0;
   double weight = tWeight;
   
   for(int i=1;i<=7; i++ ){
      sum = 0.044*pow(weight,0.75) + 1.549*expGain + sum;
      weight =  weight + expGain;  
   }      
   return(rationCost[iRt]*sum); */
      
   //Cost based on daily feed intake and estimated parameters:   
   //return( rationCost[iRt]*(dlm[iRt].k2*meanGrowth + dlm[iRt].k1[t-1]*(tWeight + 3*meanGrowth/7)  ) );
   
   //Cost based on weekly feed intake and estimated parameters:       
   double expGain = meanGrowth;
   return(rationCost[iRt]*(dlm[iRt].k1[t-1]*tWeight + dlm[iRt].k2*expGain));    
         
   //Cost based on weekly feed intake and orginal parameters:
   /* double expGain = meanGrowth;
   return(rationCost[iRt]*(7*0.044*pow(tWeight,0.75) + 1.549*expGain));*/ 
         
   //Cost based on weekly feed intake and growth parameter:
   /* return(rationCost[iRt]*2.5*expGain);*/
   
   
   //Cost based on daily feed intake and orginal parameters (version 2):       
   /*    double w=tWeight;
   if (w<=0) w=1; 
   double g=meanGrowth;
   double k1 = 1.549;
   double k2 = 0.044;
   
   return (7*g*k1 + fastPow(6*g + w,0.75)*k2 + fastPow(5*g + w,0.75)*k2 + fastPow(4*g + w,0.75)*k2 
      + fastPow(3*g + w,0.75)*k2 + fastPow(2*g + w,0.75)*k2 + fastPow(g + w,0.75)*k2 + k2*fastPow(w,0.75) ); */
         
}

// ===================================================

int HMDP::IdCountL2(int phase){
   int n, iSW, iSSd, x, t, iSG, RS;
   x=0;
   for(t = minPhaseT[phase]; t<=tMax; t++){         
   for(RS=minPhaseT[phase];RS<=t;RS++){
   if( (phase==0) & (RS!=1) ) continue;
   if(RS==tMax) continue;
   for (n=1;n<=pigs;n++){
      if ( (t<=tStartMarketing) & (n!=pigs) ) continue;
      for (iSW=0; iSW<sizeSW; iSW++){
         for (iSG=0; iSG<sizeSG; iSG++){
            for (iSSd=0; iSSd<sizeSSd; iSSd++){
               x++;
            }
         }
      }
    } 
   }
   }   
   return (x);
}

// ===================================================

int HMDP::IdCountL1(int phase){
   int n, iSW, iSSd, x, iR, RS;
   x=0;
   
   for (iR=0;iR<rations;iR++){    // rations
      if( (phase==0) & (iR!=iniFeedId) ) continue;
      for (RS=minPhaseT[phase];RS<tMax;RS++){   // possible time for changing the rations
      if((phase==0) & (RS!=1)) continue;
         for (n=1;n<=pigs;n++){   // pigs in the pen
            if( (RS<=tStartMarketing) & (n!=pigs)) continue;
            for (iSW=0; iSW<sizeSW; iSW++){    // row number in sW (C style)
               for (iSSd=0; iSSd<sizeSSd; iSSd++){  // row number in sSd
                  x++;
               }
            }
         }
      }
   }   
   return (x+2);
}

// ===================================================

void HMDP::WeightsTransPrTermTH(int & phase, int & iRation, int & t, int & n, int & iSW, int & iSG, int & iSSd) {
   if( phase!=(phases-1) ){ 
      index.assign(1, mapL1["dummy"] ); 
   }
   else{          
      index.assign(1,0); 
   }
   scope.assign(1,0); pr.assign(1,1);
   weights[0]=0; weights[1]=weightTH[t][0][iRation][iSW][iSG][iSSd][n];  // assume that iHT = 0 is the terminate threshold!
}
