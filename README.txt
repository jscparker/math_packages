
The repository contains:

   A collection of basic math routines (SVD, Eigendecomposition, FFT,
   Runge-Kutta, Random Number Generators, Arbitrary Precision Floats, etc).
   All of the routines are written in Ada. Most of the packages come
   with introductory Demo routines.

License: the software is distributed under the

   The Internet Software Consortium (ISC) License. 

   The ISC license is a permissive free software license published by the
   Internet Systems Consortium (ISC). It is functionally equivalent to the
   simplified BSD and MIT licenses, differing in its removal of language
   deemed unnecessary following the global adoption of the Berne Convention.

   The text of the license is found in file license.txt.

Optimization: a decent optimization on gcc/GNAT is provided by

   gnatmake -gnatnp -O3 -march=native -funroll-loops
     (-ffast-math is essential for best performance in golub_svd.adb)
     (-gnatNp sometimes hits a bug, evidently due to a misplaced 'pure' in GNAT)

It's a good idea to do a preliminary run which exercises Assertions,
and other Checks:

   gnatmake -Wall -gnatwa -gnatVa -gnata -gnato -fstack-check -gnateE xxx.adb

Segmentation fault: giant arrays sometimes cause segmentation fault if 
   insufficient stack mem. In Bash shell, you can try: ulimit -s unlimited.
   The -s is for stack. In csh try: limit stacksize unlimited.

directory Arbitrary contains:

   package Extended_Real:
      An arbitrary precision floating-point data type: e_Real. Lower limit
      is 28 decimals, no upper limit is enforced. The package is Pure. 
      Floating point attributes (Ada 95) are implemented as function calls.
      The package exports standard floating point operators:
      "*", "+", "/", "**", "Abs", "<", ">", "<=" , ">=", etc. 
   package Extended_Real.Elementary_Functions:
      Sin, Cos, Sqrt, Arcsin, Arccos, Arctan, Log, Exp, Reciprocal (1/x),
      Reciprocal_Nth_Root (x to the power of -1/N), Divide, and "**" for
      extended arguments and exponents. Routines are Ada 95'ish (except arctan).
   package Extended_Real.IO:
      Text to (extended precision) e_Real translation routines, and
      e_Real to Text translation routines.
   package e_Derivs:
      Extended precision routines for taking high order derivatives of
      functions.  Functions constructed from  "*", "+", "/", "**", Sin,
      Cos, Sqrt, Arcsin, Arccos, Arctan, Log, Exp, Compose = f(g(x)),
      and Reciprocal can be differentiated to order specified by user.
  
directory Disorderly contains:

   package Disorderly.Random:
      A long-period high-quality non-linear Random Number Generator.
      Procedure Get_Random produces 61 random bits per call. The package
      is pure, and run-time per call constant - designed to be suitable
      for multi-tasking sumulations.
   package Disorderly.Basic_Rand:
      a long-period high-quality linear Random Number Generator.
   package Disorderly.Random.Deviates:
      Random deviates (variates) with the following distributions:
      Uniform, Normal (Gaussian), Exponential, Lorentzian (Cauchy),
      Poissonian, Binomial, Negative Binomial, Weibull, Rayleigh, 
      Student_t, Beta, Gamma, Chi_Squared, Log_Normal, Multivariate_Normal.
   package Gamma:
      function Log_Gamma
      Natural logarithm of Gamma function for positive real arguments.
   package Chi_Gaussian_CDF:
      function Incomplete_Gamma
      function Chi_Squared_CDF
      function Normal_CDF
      The packages provides cummulative distribution functions (CDF) 
      for the Chi squared distribution, and for the standard normal 
      distribution, so that p-values can be reliably calculated
      to machine precision.
 
directory Fourier contains:

   Two routines for calculating Discrete Fourier Transforms (DFT).
   Both do the DFT in O(N*Log(N)) time.

   package Fourier8:
      Standard FFT (Radix 8); appropriate for time-sensitive tasks.
   package Chirped:
      Chirped FFT; usually best choice for data analysis.


directory linear_algebra contains:

   Several lin alg routines for real-valued matrices:

      QR_Symmetric_Eigen (Eigendecomposition for symmetric matrices.)
      Peters_Eigen       (Eigendecomposition for general square matrices.)
      Golub_SVD          (SVD decomposition for rectangular matrices.)
      Givens_QR          (QR decomposition for rectangular matrices.)
      Jacobi_Eigen       (Jacobi eigen-decomposition for symmetric matrices.)
      Cholesky_LU        (Cholesky LU decomposition for positive definite matrices.)
      Crout_LU           (Crout LU decomposition for square matrices.)
      Banded_LU          (Crout LU decomposition for banded matrices.)
      Hessenberg         (Hessenberg decomposition using Givens rotations.)
      Tridiagonal        (Tridiagonalization of symmetric matrices.)

   Several of the routines are optimized for numerical accuracy, none
   for speed. Enhanced Givens rotations are used everywhere except LU
   decompositions. Enhanced Givens rotations are sluggish, but improve
   accuracy and reliability. Reliability is generally favored over speed.

directory Ordinary contains:

   Several routines for integrating ordinary differential equations.
   packages Runge_8th and Runge_5th:
      The Prince-Dormand 8th order and 5th order Runge-Kutta integrators
      for general ordinary differential equations: (d/dt) Y(t) = F(t, Y).
      The Prince-Dormand methods offer optional error control with variable
      integration step-sizes.
   package Predictor_1:
      A 17th order Predictor-Corrector integrator for general ordinary
      differential equations: (d/dt) Y(t) = F(t, Y).
   package Predictor_2:
      A 20th order Predictor-Corrector integrator for conservative
      differential equations: (d/dt)^2 Y(t) = F(t, Y).
   Several Demo/Test routines that demonstrate the relative merits
      of the different methods.

directory polynomial/spline contains:

   package Spline:
      Routines for generating cubic spline fits to arbitrarily spaced
      data points, plus routines for interpolation, integration,
      and differentiation of the spline.
   package Tridiagonal_LU:
      Crout's method for LU decomposition of tri-diagonal matrices and
      linear equation solving. Tridiagonal_LU supports package Spline.

directory polynomial/clenshaw contains:

   package Clenshaw:
      Generates, evaluates, and sums functions defined by recurrance
      relations. Recurrance relations for generating Laguerre,
      Generalized Laguerre, Legendre, Chebychev (1st and 2nd kind),
      Associated Legendre, and Hermite polynomials are given in the
      text of package Clenshaw.
   package Chebychev:
      Uses Clenshaw's method to generate Chebychev polynomials of the
      2nd kind.
   package Spherical_Harmonics:
      Uses Clenshaw's method to generate Spherical Harmonics
      (Associated Legendre polynomials).
   package Factorial:
      Uses the Stieltjes' continued fraction method to generate natural
      logarithms of Factorials. Procedure Log_Factorial is used to
      normalize functions generated by the recurrance method.
   package Gauss_Quadrature_61:
      61 point Gauss-Kronrod quadrature rules. Gauss-Kronrod quadrature
      is used by test procedures Spherical_Harm_tst_1, Chebychev_tst_1.
   package Chebychev_Quadrature:
      Gauss-Chebychev quadrature. Gauss_Quadrature is used by test
      procedures Chebychev_tst_1 and Cheby_Quad_tst_1.

directory polynomial/least_sq_poly contains:

   package Orthogonal_Polys:
      Generates Gram-Schmidt orthogonal polynomials on arbitrarily
      spaced discrete grid points. The package also provides a routine
      for high order polynomial least-squares curve fitting.
   procedures Opolys_tst_1 and Opolys_tst_2:
      Demonstrate high order least squares polynomial fitting.

