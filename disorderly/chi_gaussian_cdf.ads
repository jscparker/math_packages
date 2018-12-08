
-------------------------------------------------------------------------------
-- package Chi_Gaussian_CDF, CDF of Normal and Chi-square distributions
-- Copyright (C) 1995-2018 Jonathan S. Parker
--
-- Permission to use, copy, modify, and/or distribute this software for any
-- purpose with or without fee is hereby granted, provided that the above
-- copyright notice and this permission notice appear in all copies.
-- THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
-- WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
-- MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
-- ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
-- WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
-- ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
-- OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
-------------------------------------------------------------------------------

-- package Chi_Gaussian_CDF
--
-- The Incomplete Gamma function for calculating area under the chi 
-- squared distribution, and also the area under the normalized
-- Gaussian with unit width (the normal distribution).
--
-- The purpose of this package is to provide the cummulative distribution 
-- functions (CDF) for the Chi squared distribution, and for the standard 
-- normal distribution, so that p-values can be reliably calculated.
-- The area under the chi squared distribution is calculated by
-- function Chi_Squared_CDF.  The area under the normal distribution
-- is calculated by function Normal_CDF.
--
--
-- Function Normal_CDF:
--
-- Evaluates the tail area of the standardized normal distribution
-- (normalized gaussian with unit width) from -inf to x.
--
-- The normal distribution:
--
--      Phi(x) = Exp (-(x-mu)^2/(2*Sigma^2))  /  Sqrt(2*Pi*Sigma^2)
--
-- has variance = Sigma^2, and mean = mu.  
-- The Standard Normal Distribution has 
-- variance = Sigma^2 = 1, and mean = mu = 0.  Phi(x) is normalized.
--
-- The cumulative normal distribution function (area under Phi(x) from -inf to x)
--
--      Phi_CDF(x) = Phi_CDF_stnd ((x-mu)/Sigma)
--
-- where Phi_CDF_stnd is the cumulative distribution function of the standard
-- normal distribution (mu=0, Sigma=1).  Notice as x -> inf, Phi_CDF -> 1 since
-- Phi is normalized.
--
-- Notes on Algorithm:
--
-- Sums a series at small argument X (see Abramowitz, Stegun 6.5.29),
-- and uses Gauss's continued fraction for large argument X.
-- Seems to be good to about 14 significant figures, when Real'Digits = 15.
--
-- Remember that this package provides a continuous distribution. If 
-- the p-values are obtained by testing discrete distributions then they
-- won't behave exactly as predicted for continuous distributions. If there
-- are a large number of degrees of freedom, you might need (for example) to 
-- divide interval over which statistic is measured into (say) 2**10 to 2**20 
-- segments to get discrete behaviour to approach that of continuous
-- distribution.
--
--
--  Inc_Gamma(a, w) = Int{from 0 to w} [exp(-t) t**(a-1) dt]  / Gamma(a)
--
--  Gamma(0.5) = Sqrt(pi) and let t = u*u/2:
--
--  Inc_Gamma(0.5, w) = Int{from 0 to sqrt(2w)} [2 exp(-v*v/2) dv]  / Sqrt(2 pi)
--
--  Inc_Gamma(0.5, x*x/2) / 2 = Int{from 0 to x)} [exp(-v*v/2) dv]  / Sqrt(2 pi)
--

generic

  type Real is digits <>;
  -- Anything up to 16 digit floats is appropriate here; 18 is probably OK.
  -- Haven't tried 32 digit floats, tho they should work in principal.

package Chi_Gaussian_CDF is

  function Incomplete_Gamma_C (a : Real; x : Real) return Real;

  function Incomplete_Gamma (a : Real;  x : Real)  return Real;

  function Chi_Squared_CDF 
    (True_Degrees_Freedom : Real; 
     Chi_Squared          : Real) 
     return Real;
  --  Calculate the Chi-Squared statistic, plug it into this function,
  --  and you get a p-value, a distribution that is uniformly distributed
  --  on [0, 1] if the statistic follows the standard Chi-squared probability
  --  distribution.

  function Normal_CDF (x : Real) return Real;

  procedure Test_Normal_CDF (x_0 : in Real; Error : out Real);

end Chi_Gaussian_CDF;
