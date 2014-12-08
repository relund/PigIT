#include "cumNorm.h"




//' Univariate cumulative normal (distribution function)
//' 
//' Based on the paper "Better approximations to cumulative normal functions" by 
//' G. West \url{http://www.wilmott.com/pdfs/090721_west.pdf}. The source code
//' has been taken from the authors webpage \url{http://finmod.co.za/research.html}.
//' 
//' @param x Value of quantiles.
//' @param mean The mean.
//' @param sd The standard deviation.
//' @return The distribution function.
//' @author Lars Relund \email{lars@@relund.dk}
//' @examples
//' # test against pnorm
//' times<-500
//' x<-runif(times,-10,10)
//' mean<-0
//' sd<-3
//' diff<-abs(pnorm(x,mean,sd)-pNorm1D(x,mean,sd))
//' if (all(diff<1e-15)) cat("No problems")
//' 
//' # compare performance
//' testpnorm<-function() {
//'    x<-runif(1,-10,10)
//'    mean<-runif(1,-10,10)
//'    sd<-runif(1,0,5)
//'    pnorm(x,mean,sd)
//' }
//' testpNorm1D<-function() {
//'    x<-runif(1,-10,10)
//'    mean<-runif(1,-10,10)
//'    sd<-runif(1,0,5)
//'    pNorm1D(x,mean,sd)
//' }
//' 
//' require(rbenchmark)
//' benchmark(testpnorm(), testpNorm1D(), order="relative", replications=100000)[,1:4]
//' @export
// [[Rcpp::export]]
SEXP pNorm1D(NumericVector x, double mean, double sd) {
   int n = x.size();
   NumericVector out(n);
   for(int i = 0; i < n; ++i) {
      out[i] = cumNorm1D( (x[i]-mean)/sd );
   }
   return out;
}








//' Bivariate cumulative normal (distribution function) using Armadillo
//' 
//' Based on the paper "Better approximations to cumulative normal functions" by 
//' G. West \url{http://www.wilmott.com/pdfs/090721_west.pdf}. The source code
//' has been taken from the authors webpage \url{http://finmod.co.za/research.html}.
//' 
//' Transform into a bivariate with sd = 1 using Y = 1/sqrt(diag(sigma))(X-mean), i.e.
//' Pr(X<=x) = P(Y<=1/sqrt(diag(sigma))(x-mean)).
//' 
//' @param lower The vector of lower limits (2-dim vector).
//' @param upper The vector of upper limits (2-dim vector).
//' @param mean The mean vector (2-dim vector).
//' @param sigma The covariance matrix (2x2). 
//' @return The distribution function.
//' @author Lars Relund \email{lars@@relund.dk}
//' @examples
//' # test against pmvnorm
//' times<-500
//' diff<-rep(NA,times)
//' diff_arma<-rep(NA,times)
//' for (i in 1:times) {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    p1<-pmvnorm(lower=lower,upper=upper,mean=mean,sigma=sigma)
//'    p2<-pNorm2D(lower, upper, mean, sigma)
//'    p3<-pNorm2D_arma(lower, upper, mean, sigma)
//'    diff[i]<-abs(p1-p2)
//'    diff_arma[i]<-abs(p1-p3)
//' }
//' if (all(diff<1e-15)) cat("No problems with pNorm2D") else stop("Results for pNorm2D not the same!")
//' if (all(diff_arma<1e-15)) cat("No problems with pNorm2D_arma") else stop("Results for pNorm2D_arma not the same!")
//' 
//' # compare performance
//' testPmvnorm<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pmvnorm(lower=lower,upper=upper,mean=mean,sigma=sigma)
//' }
//' testpNorm2D<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pNorm2D_arma(lower, upper, mean, sigma)
//' }
//' testpNorm2D_arma<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pNorm2D_arma(lower, upper, mean, sigma)
//' }
//' 
//' require(rbenchmark) 
//' benchmark(testPmvnorm(), testpNorm2D(), testpNorm2D_arma, order="relative", replications=10000)[,1:4]
//' @export
// [[Rcpp::export]]
double pNorm2D_arma(arma::vec lower, arma::vec upper, arma::vec mean, arma::mat sigma) {
   arma::vec lb = sigma.diag();
   arma::vec ub = sigma.diag();
   double rho;
   
   lb = (lower-mean)/sqrt(lb);
   ub = (upper-mean)/sqrt(ub);
   rho = sigma(1,0)/sqrt(sigma(0,0)*sigma(1,1));
   //Rcout << lb << " " << ub << " " << rho << endl;
   
   double p1 = cumNorm2D(ub[0], ub[1], rho);
   double p2 = cumNorm2D(lb[0], lb[1], rho);
   double p3 = cumNorm2D(ub[0], lb[1], rho);
   double p4 = cumNorm2D(lb[0], ub[1], rho);
   //Rcout << p1 << " " << p2 << " " << p3 << " " << p4 << endl;
   return fmax(0,p1-p3-p4+p2);   // use max since sometimes may be negative
}



//' Bivariate cumulative normal (distribution function)
//' 
//' Based on the paper "Better approximations to cumulative normal functions" by 
//' G. West \url{http://www.wilmott.com/pdfs/090721_west.pdf}. The source code
//' has been taken from the authors webpage \url{http://finmod.co.za/research.html}.
//' 
//' Transform into a bivariate with sd = 1 using Y = 1/sqrt(diag(sigma))(X-mean), i.e.
//' Pr(X<=x) = P(Y<=1/sqrt(diag(sigma))(x-mean)).
//' 
//' @param lower The vector of lower limits (2-dim vector).
//' @param upper The vector of upper limits (2-dim vector).
//' @param mean The mean vector (2-dim vector).
//' @param sigma The covariance matrix (2x2). 
//' @return The distribution function.
//' @author Lars Relund \email{lars@@relund.dk}
//' @examples
//' # test against pmvnorm
//' times<-500
//' diff<-rep(NA,times)
//' diff_arma<-rep(NA,times)
//' for (i in 1:times) {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    p1<-pmvnorm(lower=lower,upper=upper,mean=mean,sigma=sigma)
//'    p2<-pNorm2D(lower, upper, mean, sigma)
//'    p3<-pNorm2D_arma(lower, upper, mean, sigma)
//'    diff[i]<-abs(p1-p2)
//'    diff_arma[i]<-abs(p1-p3)
//' }
//' if (all(diff<1e-15)) cat("No problems with pNorm2D") else stop("Results for pNorm2D not the same!")
//' if (all(diff_arma<1e-15)) cat("No problems with pNorm2D_arma") else stop("Results for pNorm2D_arma not the same!")
//' 
//' # compare performance
//' testPmvnorm<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pmvnorm(lower=lower,upper=upper,mean=mean,sigma=sigma)
//' }
//' testpNorm2D<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pNorm2D_arma(lower, upper, mean, sigma)
//' }
//' testpNorm2D_arma<-function() {
//'    sigma <- matrix(c(runif(1,1,10),1,1,runif(1,1,10)),nrow=2)
//'    sigma[1,2]<-sigma[2,1]<-sqrt(sigma[1,1]*sigma[2,2])*runif(1,-1,1)
//'    mean <- c(runif(1,1,10),runif(1,1,10))
//'    lower<-c(runif(1,-6,6),runif(1,-6,6))
//'    upper<-c(runif(1,lower[1],6),runif(1,lower[2],6))
//'    pNorm2D_arma(lower, upper, mean, sigma)
//' }
//' 
//' require(rbenchmark) 
//' benchmark(testPmvnorm(), testpNorm2D(), testpNorm2D_arma, order="relative", replications=10000)[,1:4]
//' @export
// [[Rcpp::export]]
double pNorm2D(NumericVector lower, NumericVector upper, NumericVector mean, NumericMatrix sigma) {
   NumericVector lb(2);
   NumericVector ub(2);
   double rho;

   for (int i=0;i<2;i++) {
      lb[i] = (lower[i] - mean[i])/sqrt(sigma(i,i));
      ub[i] = (upper[i] - mean[i])/sqrt(sigma(i,i));
   }
   rho = sigma(1,0)/sqrt(sigma(0,0)*sigma(1,1));
   
   double p1 = cumNorm2D(ub[0], ub[1], rho);
   double p2 = cumNorm2D(lb[0], lb[1], rho);
   double p3 = cumNorm2D(ub[0], lb[1], rho);
   double p4 = cumNorm2D(lb[0], ub[1], rho);
   //Rcout << p1 << " " << p2 << " " << p3 << " " << p4 << endl;
   return fmax(0,p1-p3-p4+p2);   // use max since sometimes may be negative
}

