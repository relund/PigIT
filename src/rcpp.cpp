#include "hmdp.h"

using namespace Rcpp;
using namespace std;

//' Build the HMDP (3 levels with shared linking) using the C++ binary writer. 
//' 
//' The MDP specified at level 3 is generated for each ration for a dummy state at the first stage at level 2. 
//' 
//' @param filePrefix Prefix used by the binary files storing the MDP model.
//' @param param Model parameters a list created using \code{\link{setParameters}}.
//' @param dlms A list a long as the number of rations where each entry is a DLM created using \code{\link{iniDLM}} and 
//'   \code{\link{buildDLM}}.
//' @param dglmParam Model parameters a list created using \code{\link{iniDGLM}}.
//'
//' @return Build log (character).
//' @author Lars Relund \email{lars@@relund.dk}
//' @export
// [[Rcpp::export]]
SEXP BuildHMDP(const CharacterVector filePrefix, const List param, const List dlms, const List dglmParam) {
   string prefix = as<string>(filePrefix);
   HMDP Model(prefix, param, dlms, dglmParam);
   Rcout << "Total number of states: " << Model.countStatesHMDP() << endl;
   return( Model.BuildHMDP() );
   //return(wrap(0));
}















