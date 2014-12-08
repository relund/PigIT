#ifndef CUMNORM_H
#define CUMNORM_H

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <cmath>
using namespace Rcpp;
using namespace std;


/*
' Univariate cumulative normal (distribution function for the standard normal)
' 
' Based on the paper "Better approximations to cumulative normal functions" by 
' G. West \url{http://www.wilmott.com/pdfs/090721_west.pdf}. The source code
' has been taken from the authors webpage \url{http://finmod.co.za/research.html}.
' 
' @param x Value of quantile (number).
' @return The distribution function.
' @author Lars Relund \email{lars@@relund.dk}
*/
inline double cumNorm1D(double x)
{
   double cn = 0;
	
	double xabs = fabs(x);
	if (xabs > 37) 
	{
		cn = 0;
	}
	else
	{	double exponential = exp(-pow(xabs,2)/2);
		if (xabs < 10 * M_SQRT1_2) 
		{
			double build = 3.52624965998911E-02 * xabs + 0.700383064443688;
			build = build * xabs + 6.37396220353165;
			build = build * xabs + 33.912866078383;
			build = build * xabs + 112.079291497871;
			build = build * xabs + 221.213596169931;
			build = build * xabs + 220.206867912376;
			cn = exponential * build;
			build = 8.83883476483184E-02 * xabs + 1.75566716318264;
			build = build * xabs + 16.064177579207;
			build = build * xabs + 86.7807322029461;
			build = build * xabs + 296.564248779674;
			build = build * xabs + 637.333633378831;
			build = build * xabs + 793.826512519948;
			build = build * xabs + 440.413735824752;
			cn = cn / build;
		}
		else
		{
			double build = xabs + 0.65;
			build = xabs + 4.0 / build;
			build = xabs + 3.0 / build;
			build = xabs + 2.0 / build;
			build = xabs + 1.0 / build;
			cn = exponential / build / 2.506628274631;
		}
	}

	if (x > 0) 
	{
		cn = 1 - cn;
	}
	return cn;
}





// constants used for bivariate calculation
const double XX[10][3] =
{
   {-0.932469514203152,-0.981560634246719,-0.993128599185095},
	{-0.661209386466265,-0.904117256370475,-0.963971927277914},
	{-0.238619186083197,-0.769902674194305,-0.912234428251326},
	{0.,-0.587317954286617,-0.839116971822219},
	{0.,-0.36783149899818,-0.746331906460151},
	{0.,-0.125233408511469,-0.636053680726515},
	{0.,0.,-0.510867001950827},
	{0.,0.,-0.37370608871542},
	{0.,0.,-0.227785851141645},
	{0.,0.,-0.0765265211334973}
};

const double W[10][3] =
{
	{0.17132449237917,0.0471753363865118,0.0176140071391521},
	{0.360761573048138,0.106939325995318,0.0406014298003869},
	{0.46791393457269,0.160078328543346,0.0626720483341091},
	{0.,0.203167426723066,0.0832767415767048},
	{0.,0.233492536538355,0.10193011981724},
	{0.,0.249147045813403,0.118194531961518},
	{0.,0.,0.131688638449177},
	{0.,0.,0.142096109318382},
	{0.,0.,0.149172986472604},
	{0.,0.,0.152753387130726}
};

/*
' Bivariate cumulative normal (distribution function for the standard normal (using a correlation matrix))
' 
' Based on the paper "Better approximations to cumulative normal functions" by 
' G. West \url{http://www.wilmott.com/pdfs/090721_west.pdf}. The source code
' has been taken from the authors webpage \url{http://finmod.co.za/research.html}.
' 
' @param x Value of first quantile (double).
' @param y Value of second quantile (double).
' @param correlation cofficient.
' @return The distribution function.
' @author Lars Relund \email{lars@@relund.dk}
*/
inline double cumNorm2D(double x, double y, double correlation)
{
	int NG;
	int LG;
	
	if (fabs(correlation) < 0.3)
	{
		NG = 1;
		LG = 3;
	}
	else if (fabs(correlation) < 0.75)
	{
		NG = 2;
		LG = 6;
	}
	else 
	{
		NG = 3;
		LG = 10;
	}

	double h = -x;
	double k = -y;
	double hk = h * k;
	double BVN = 0;

	if (fabs(correlation) < 0.925)
	{
		if (fabs(correlation) > 0)
		{
		    double hs = (h * h + k * k) / 2;
		    double asr = asin(correlation);
			for (int i = 1; i <= LG; ++i)
			{
				for (int iss = -1; iss <=1; iss += 2)
				{
					double sn = sin(asr * (iss * XX[i-1][NG-1] + 1) * 0.5);
					BVN = BVN + W[i-1][NG-1] * exp((sn * hk - hs) / (1.0 - sn * sn));
				}
			}
			BVN = BVN * asr * 0.795774715459476678e-1;
		}
		BVN = BVN + cumNorm1D(-h) * cumNorm1D(-k);
	}
	else
	{
		if (correlation < 0) 
		{
			k *= -1;
			hk *= -1;
		}
		if (fabs(correlation) < 1)
		{
		    double Ass = (1 - correlation) * (1 + correlation);
			double a = sqrt(Ass);
			double bs = (h-k)*(h-k);
			double c = (4 - hk) / 8;
			double d = (12 - hk) / 16;
			double asr = -(bs / Ass + hk) / 2;
			if (asr > -100)
			{
				BVN = a * exp(asr) * (1 - c * (bs - Ass) * (1 - d * bs / 5) / 3 + c * d * Ass * Ass / 5);
			}
			if (-hk < 100)
			{
				double B = sqrt(bs);
				BVN = BVN - exp(-hk / 2) * 2.506628274631 * cumNorm1D(-B / a) * B * (1 - c * bs * (1 - d * bs / 5) / 3);
			}
			a /= 2;
			for (int i = 1; i <= LG; ++i)
			{
				for (int iss = -1; iss <= 1; iss += 2)
				{
					double xs = a * (iss * XX[i-1][NG-1] + 1);	
					xs = fabs(xs*xs);
					double rs = sqrt(1 - xs);
					asr = -(bs / xs + hk) / 2;
					if (asr > -100)
					{
						BVN = BVN + a * W[i-1][NG-1] * exp(asr) * (exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs - (1 + c * xs * (1 + d * xs)));
					}	
				}
			}
			BVN *= - 0.159154943091895336;
		}
		if (correlation > 0)
		{
			BVN = BVN + cumNorm1D(-max(h, k));
		}
		else
		{
			BVN *= -1;
			if (k > h)
			{
				BVN = BVN + cumNorm1D(k) - cumNorm1D(h);
			}
		}
	}	 
	return BVN;
}

double pNorm2D_arma(arma::vec lower, arma::vec upper, arma::vec mean, arma::mat sigma);


#endif




