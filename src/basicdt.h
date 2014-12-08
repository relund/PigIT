#ifndef BASICDT_HPP
#define BASICDT_HPP

//#include <iostream>
//#include <string>
#include <sstream>
#include <vector>
using namespace std;

// Basic datatypes and definitions
//-----------------------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )  ///< If x smaller than y then return x, otherwise return y
#define MAX(x,y) ( (x) > (y) ? (x) : (y) )  ///< If x larger than y then return x, otherwise return y

typedef double flt;                 ///< A floating number datatype.
typedef double* fltPtr;             ///< A floating number pointer.
typedef unsigned int idx;           ///< A datatype for storing array indexes etc. Note don't make it short since also used as index for states! //unsigned since then have to cast!!
typedef idx* idxPtr;                ///< Pointer to idx.
typedef unsigned short int uSInt;   ///< Unsigned short integer.
typedef unsigned int uInt;          ///< Unsigned long integer
typedef int* intPtr;


const flt INF=18000000000000000.0;      ///< Infinity (or values above).
const uInt INFINT = 1000000000;         ///< Infinity integer.
const flt PRECISION = 1e-10;  ///< used for comparison floats
//const int FALSE = 0;
//const int TRUE = 1;

/** Global function for comparing two floats. Assume equal if their difference
* if less than PRECISION.
* \return True if equal.
*/
inline bool Equal(flt n1,flt n2) {
    return ((n2-PRECISION)<=n1 && n1<=(n2+PRECISION));
};

/** Global function for comparing two floats. Assume equal if their difference
* if less than precision.
* \return True if equal.
*/
inline bool Equal(flt n1, flt n2, flt precision) {
    return ((n2-precision)<=n1 && n1<=(n2+precision));
};

inline bool Equal(int n1,int n2) {
    return (n1==n2);
};

/** Global function for converting a number to a string */
template <typename T>
std::string inline ToString(T t) {
    std::ostringstream s;
    s << t;
    return s.str();
};

/** Global function for converting a string
 \param t The variable of the result.
 \param s The string.
 \param *f One of std::hex, std::dec or std::oct.
 \return 0 if failed and 1 otherwise.
 Example

 if(from_string<float>(f, std::string("123.456"), std::dec))
  {
    std::cout << f << std::endl;
  }
  else
  {
    std::cout << "from_string failed" << std::endl;
  }
 */
template <typename T>
bool inline from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
};


/** convert a vector to a comma separated string. */
template <typename T>
string inline vec2String(const vector<T>& v) {
    std::ostringstream s;
    s << "(";
    for (idx i=0; i<v.size()-1; ++i) s << v[i] << ",";
    s << v[v.size()-1] << ")";
    return s.str();
}



/** Global function for converting a flt to a string */
/*inline string flt2String(const flt i)
{
  ostringstream stream;
  stream << i;
  return stream.str();
};*/
//-----------------------------------------------------------------------------
//#define null 0
//-----------------------------------------------------------------------------
  
/** Global function for approximating the area under the curve  
 * @param totalLength The length of the considered interval. 
 * @param y A vector for the average values of the heights of the divided rectangular.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline Integral(const double & totalLength, const arma::vec & y) {
   return ( totalLength/y.size() ) * sum(y);
}


/** Global function for computing the one side trucated expected velue (upper tail in normal distribution) 
  * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline ExpUTail(const double & mean, const double & sigma, const double & uTail) {
   double alpha=(uTail-mean)/sigma;
   double lambda=(double) (R::dnorm(alpha,0,1,0))/( (double) (1-R::pnorm(alpha,0,1,1,0))  ); 
   return(mean+sigma*lambda);
}


/** Global function for computing the one side trucated variance(upper tail in normal distribution)
  * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline VarUTail(const double & mean, const double & sigma, const double & uTail) {
   double alpha=(uTail-mean)/sigma;
   double lambda=(double) (R::dnorm(alpha,0,1,0))/( (double) (1-R::pnorm(alpha,0,1,1,0))  ); 
   double Xi=lambda*(lambda-alpha);
   return(pow(sigma,2)*(1-Xi));
}


/** Global function for computing the one side trucated expected velue (lower tail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline ExpLTail(const double & mean, const double & sigma, const double & lTail) {
   double beta=(lTail-mean)/sigma;
   double coe= (double) (R::dnorm(beta,0,1,0))/( (double) (R::pnorm(beta,0,1,1,0))  ); 
   return(mean - sigma*coe);
}


/** Global function for computing the one side trucated variance (lower tail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline VarLTail(const double & mean, const double & sigma, const double & lTail) {
   double beta=(lTail-mean)/sigma;
   double coe= (double) (R::dnorm(beta,0,1,0) )/( (R::pnorm(beta,0,1,1,0))  ); 
   return(pow(sigma,2)*(1-beta*coe-pow(coe,2)));
}


/** Global function for computing two sided trucated expected velue (in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline ExpTruncated(const double & mean, const double & sigma, const double & uTail, const double & lTail) {
   double alpha=(uTail-mean)/sigma;
   double beta=(lTail-mean)/sigma;
   double coe= (double) (R::dnorm(alpha,0,1,0) - R::dnorm(beta,0,1,0))/( (double) (R::pnorm(beta,0,1,1,0) - R::pnorm(alpha,0,1,1,0))   );
   return(mean + sigma*coe) ;
}


/** Global function for computing two sides trucated variance (in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline VarTruncated(const double & mean, const double & sigma, const double & uTail, const double & lTail) {
   double alpha=(uTail-mean)/sigma;
   double beta=(lTail-mean)/sigma;
   double coe= (double) (R::dnorm(alpha,0,1,0) - R::dnorm(beta,0,1,0))/( (double) (R::pnorm(beta,0,1,1,0) - R::pnorm(alpha,0,1,1,0)) );
   double coe1= (double) (alpha*R::dnorm(alpha,0,1,0) - beta*R::dnorm(beta,0,1,0))/( (double) (R::pnorm(beta,0,1,1,0) - R::pnorm(alpha,0,1,1,0)) );
   return(pow(sigma,2)*(1+coe1-pow(coe,2))) ;
}


/** Global function for computing two sides trucated of E(X^2) (in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline ExpX2Truncated(const double & mean, const double & sigma, const double & uTail, const double & lTail) {
   
   return( VarTruncated(mean, sigma, uTail, lTail) + pow(ExpTruncated(mean, sigma, uTail, lTail),2) );

}


/** Global function for computing one sides trucated of E(X^2) (Utail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline ExpX2UTail(const double & mean, const double & sigma, const double & uTail) {
   
   return( VarUTail(mean, sigma, uTail) + pow(ExpUTail(mean, sigma, uTail),2) );

}


/** Global function for computing the distribution function for Utail (F(x|x>Utail) Utail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline cdfuTail(double x, const double & mean, const double & sigma, const double & uTail) {
   
   double alpha=(uTail-mean)/sigma;
   double delta=(x-mean)/sigma;
   double Z=1-R::pnorm(alpha,0,1,1,0);
   return((double) (R::pnorm(delta,0,1,1,0) - R::pnorm(alpha,0,1,1,0))/(Z));
}


/** Global function for computing the probability distribution function for Ltail (P(x|x<lTail) Utail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline pdflTail(double x, const double & mean, const double & sigma, const double & lTail) {
   
   double beta=(lTail-mean)/sigma;
   double delta=(x-mean)/sigma;
   double Z=R::pnorm(beta,0,1,1,0);
   return( (double) (R::dnorm(delta,0,1,0) ) /(Z*sigma) );
}

/** Global function for computing the distribution function for Ltail (F(x|x<lTail) Utail in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline cdflTail(double x, const double & mean, const double & sigma, const double & lTail) {
   
   double prob;
   double beta=(lTail-mean)/sigma;
   double delta=(x-mean)/sigma;
   double Z=R::pnorm(beta,0,1,1,0);
   if(lTail<x){
     prob=1; 
   } 
   else{
      prob=(double) (R::pnorm(delta,0,1,1,0) )/(Z);
   }
   return( prob );
}


/** Global function for computing the distribution function for  F(x|uTail<x<lTail)  in normal distribution) 
 * @param mean Mean of normal distribution.
 * @param sigma Standard deviation of the normal distribution.
 * @param lTail Lower tail of the truncated normal distribution.
 * @param uTail Upper tail of the truncated normal distribution.
 * 
 * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
 */
double inline cdfTail(double x, const double & mean, const double & sigma, const double & uTail, const double & lTail) {
   
   double alpha=(uTail-mean)/sigma;
   double beta=(lTail-mean)/sigma;
   double delta=(x-mean)/sigma;
   double Z=R::pnorm(beta,0,1,1,0) - R::pnorm(alpha,0,1,1,0);
   return( (double) (R::pnorm(delta,0,1,1,0) - R::pnorm(alpha,0,1,1,0)) /(Z) );
}


/** Calculate a^b
 * @param a Base.
 * @param b Exponent
 * 
 * @author See http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
 */
double inline fastPow( const double & a,  const double & b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}


#endif
