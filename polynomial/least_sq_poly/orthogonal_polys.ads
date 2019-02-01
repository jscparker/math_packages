
-- package Orthogonal_Polys
--
-- Generates orthogonal polynomials on dicrete grid points using the
-- Gram-Schmidt method. The package also provides a routine for
-- polynomial least-squares curve fitting: procedure Poly_Fit. 
--
-- The least-squares curve fitting routine could have been defined in
-- another package, but it turned out to be useful in the test procedures,
-- and least-squares curve fitting is the primary application of the code.
-- The algorithm is due to Shampine (1), who was the first to show how
-- high order least squares polynomial fits could be performed on arbitrary
-- data sets. 
--
-- Given a set of X-axis grid points {X_0, ..., X_m} and a set of weights
-- associated with those points {W_0, ..., W_m} there is a unique set of
-- Gram-Schmidt orthogonal polynomials {Q_0(X), ..., Q_m(X)} that range in
-- degree from 0 to m, where m is one less than the number of distinct
-- data points.
--
-- This package defines routines for generating the polynomials Q_j(X)
-- and routines for calculating values, derivatives, and integrals of
-- functions F(X) composed of linear combinations of these polynomials:
-- F(X) = SUM [C_j * Q_j(X)].  The Clenshaw summation formula is used
-- for most of the work.
--
-- The polynomials generated are Legendre polynomials if the weights
-- are constant.  If the weights go like W_j = Sqrt (1 - X_j**2) on
-- X in [-1,1], then Tchebychev polys are generated.  Other well-known
-- polys are associated with other weights and intervals on X.  These
-- polys are not identical to the kind associated with continuous
-- intervals of X, (although they agree in the limit of large data sets.)
-- The continuous kind are orthogonal respect integration:
-- Integral_Of [Q_j(X) * Q_k(X) * W(X) * dX] = 0 if j /= k.  
-- The discrete
-- polys are orthogonal respect an inner product, a summation over
-- discrete data points SUM [Q_j(X_i) * Q_k(X_i) * W_i] = 0 if j /= k.
-- (A different package provides routines that generate values of the
-- continuous polys.  It also uses the Clenshaw summation formula.)
--
--  Notes On Use
--
--  The procedure Poly_Fit fits a polynomial in variable X to a
--  function F(X) defined by a set of points in the form (X_i, F(X_i)).
--  The fit minimizes the sum of the squares of the residuals,
--  the error terms Error_i in the equation:
--
--                 N
--      F (X_i) = SUM ( C_m * Q_m(X_i) ) + Error_i
--                m=0
--
--  where the functions Q are orthogonal polynomials, and N is the
--  order of the polynomial we are fitting to F.
--  A linear combination of Q_m, with coefficients C_m is the best fit:
--
--                          N
--   Best_Fit_Poly (X)  =  SUM ( C_m * Q_m (X) )
--                         m=0
--
--  The procedure returns this sum, and also the coefficients C_m.
--
--  C(i) is the array of coefficients of the orthogonal Polynomials
--  that are fit to the data.  Alpha(i) and Beta(i) are returned in the
--  record Poly_Set of type Polynomial.  These are
--  used to generate the orthogonal polynomials via the recurrence relations
--  given below.
--
--  If the user wants the values of the best fit polynomial at values
--  of X not among those in the data set, then one uses function
--  Poly_Value to sum the polynomial at X via the recurrence relation
--  defined by Alpha, Beta and C.
--
--  Often one wants the coefficients of the powers of X in the
--  least-squares polynomial, instead of the coefficients of the
--  orthogonal polynomials Q_m.  For example, one commonly want
--  the best linear fit to the data, or the best parabolic fit
--  to the data,  F = E_0 + E_1 * X + E_2 * X**2,
--  and so one wants the coefficients E_m in the sum:
--
--                          N
--   Best_Fit_Poly (X)  =  SUM ( E_m * X**m ).
--                         m=0
--
--  The coefficients E_m are returned in the array Poly_Coeffients by
--  the procedure Get_Coeffs_Of_Powers_Of_X.
--
--  The procedure Poly_Fit also returns:
--
--    Mean_Square_Error = SUM (Error_i * Error_i) / Number_Of_Data_Points.
--
--  If the data set is noisy, then it often happens that as you raise the
--  order polynomial, the Mean_Square_Error quickly reaches some minimum.
--  On the other hand, if you are fitting to a noiseless function, say
--  to calculate its derivatives numerically, then one can often see
--  continued improvement in the fit even as the order of the polynomial
--  becomes very large.
--
-- Notes On the Algorithm
--
-- Take any set of orthogonal functions, {Q_n}, project them one by one
-- onto a data set exactly as you would construct a Fourier series,
-- and you have performed a linear least squares fit to the data.
-- More precisely, you have found coefficients C_n such that the
-- function formed by summing the orthogonal functions Q_n with
-- coefficients C_n has minimum distance from the data in the least
-- squares sense.  No matrix inversion is required because the
-- functions Q_n are orthogonal.
-- A more rigorous definition:  Given functions
-- Q_0, Q_1, ... Q_n, ...we want to find coefficients C_0, C_1, ...
-- C_n,... such that the sum of the squares of the errors Error_i are
-- minimized in the equation:
--
--                  N
--       F (X_i) = SUM ( C_m * Q_m(X_i) ) + Error_i,
--                 m=0
--
-- where the functions Q will be orthogonal in this package, and N is the
-- order of the polynomial we are fitting to F.  The problem is best
-- formulated in terms of inner products with respect weights
-- W(X_i). Inner product is defined (A, B) = SUM (A(X_i)*B(X_i)*W(X_i)).
-- You get the formula for the least squares fit by setting the gradiant
-- of the sum of the squares of Error_i to zero and solving for the
-- coefficients C_n.  The equation to solve for C_n turns out to be:
--
--                      N
--          (Q_n, F) = SUM ( C_m * (Q_n, Q_m) ).
--                     m=0
--
-- This is a matrix equation whose solution requires an LU decomposition
-- or the like when the basis functions Q are not orthogonal.  But
-- when the basis functions Q are orthogonal, the solution of this
-- set of simultaneous equations is trivial:
--
--           C_n = (Q_n, F) / (Q_n, Q_n),
--
-- since the matrix on the RHS of the matrix equation is diagonal.
--
-- Here is the algorithm.  First we construct the orthogonal polynomials
-- Q with respect to the weights W_i as follows:
--
--          Q_0 = 1
--          Q_m = x*Q_m-1 - SUM (B_mj * Q_j )
--                          j<k
--
-- The sum on the RHS is there only to make Q_m orthogonal to all previous
-- Q's.  The coefficients Beta_mj that achieve this are
--
--         Beta_mj = (x*Q_m-1, Q_j) / (Q_j, Q_j).
--
-- To see this take the inner product of the equation for Q_m with Q_j,
-- j<k, set it equal to zero, use recursion. The form given above is the
-- one that generalizes to many dimensions.  In one dimension the formula
-- simplifies so that only two of the Beta coefficients are nonzero.
-- We call these two Alpha and Beta:
--
--                Q_0  =  1
--                Q_1  = (X - Alpha_1)
--                Q_m  = (X - Alpha_m) * Q_m-1  -  Beta_m*Q_m-2
-- where
--                Alpha_m  =  (X*Q_m-1, Q_m-1) / (Q_m-1, Q_m-1)
--                Beta_m   =  (X*Q_m-1, Q_m-2) / (Q_m-2, Q_m-2)
--
-- Can be shown:  Beta_m   =    (Q_m-1, Q_m-1) / (Q_m-2, Q_m-2) which is
-- the form used below. (Just use X*Q_m-2 = (X-Alpha)*Q_m-2 + ... = Q_m-1 + ...).
--
-- function Poly_Value:
--
-- So now we've calculated coefficients A_m and B_m which, via the
-- recurrence relation given above, generate the orthogonal polynomials
-- Q_m.  Along the way we calculated the coefficients C_m which
-- sum the Q_m's to give a least-squares fit to the function F. Now
-- to make use of these coefficients to calculate the value of the
-- the least-squares polynomial at any X, we use Clenshaw's formula.
-- Clenshaw's method is a well known algorithm for summing any
-- series C_m*Q_m(X) where Q_m is defined by a recurrence relation.
-- Clenshaw's recurrence formula is used in everything that follows,
-- so we review it next in its most general form.
--
-- Clenshaw's formula is an algorithm for evaluating sums over
-- functions Q_m (X) that are defined by recurrence relations of the sort:
--
--            Q_0  =  some given function of X
--            Q_1  =  Alpha (1, X) * Q_0 (X)
--            Q_m  =  Alpha (m, X) * Q_m-1 (X) + Beta (m, X) * Q_m-2 (X)
--
-- The procedure here evaluates the sum
--                    n
--          F_n(X) = SUM ( C_m * Q_m(X) )
--                   m=0
--
-- where the coefficients C_m are given quantities, as are Alpha and
-- Beta. The scheme implemented here is stabler and more accurate
-- than attempts to directly perform the sum for F_n.
--
-- Clenshaw proved that F_n is given by the following formula
--
--   (*)      F_n(X) = D_0*Q_0(X) + D_1*(Q_1(X) - Alpha(1,X)*Q_0(X)).
--
-- where the D_m are functions of X that satisfy:
--
--      D_n+2 = 0
--      D_n+1 = 0
--      D_m   = C_m + Alpha(m+1,X) * D_m+1(X) + Beta(m+2,X) * D_m+2(X)
--
-- The proof of (*) is straightforward. Solve for C_m in the equation for
-- D_m above and plug it into the sum that defines F_n to get
--
--           n
--    F_n = SUM (D_m - Alpha(m+1)*D_m+1 - Beta(m+1)*D_m+2 ) * Q_m
--          k=0
--                                 n
--    F_n = D_0*Q_0 + D_1*Q_1  +  SUM ( D_m*Q_m )
--                                m=2
--                                 n
--           -Alpha(1)*D_1*Q_0 +  SUM (-Alpha(m)*D_m*Q_m-1 )
--                                m=2
--                                 n
--                             +  SUM (-Beta(m) *D_m*Q_m-2 ).
--                                m=2
--
-- Now factor out D_m from the three SUM terms above, and notice
-- that what remains is just the recurrance relation that defines
-- Q_m for m > 1.  It evaluates to zero, leaving
--
--        F_n(X) = D_0*Q_0(X) + D_1*(Q_1(X) - Alpha(1,X)*Q_0(X)).
--
-- In the special case of Othogonal polynomials, where Q_0 = 1, and
-- Q_1(X) = Alpha(1,X)*Q_0(X), we have that F_n(X) = D_0.
--
--
-- procedure Poly_Derivatives:
--
-- How do we get the derivatives of the best-fit polynomial?  Just
-- take the derivative of the Clenshaw recurrence formula given above.
-- In the special case of orthogonal polynomials it is particulary
-- easy.  By differentiating the formula given above p times it's easy
-- to see that the p-th derivative of the D_m functions of X satisfy:
--
--   D_n+2(p,X) = 0
--   D_n+1(p,X) = 0
--   D_m(p,X)
--     = p*D_m+1(p-1,X) + (X - Alpha_m)*D_m+1(p,X) - Beta_m*D_m+2(p,X)
--
-- where D(p,X) is the pth derivative of D(X).  It follows that the
-- p-th derivative of the sum of C_m*Q_m(X) equals D_0(p,X).
--
--
-- procedure  Get_Coefficients_Of_Powers_Of_X:
--
-- Clenshaw's formula is used to get Coeffs.  The
-- coefficients of the following D polynomials are put into arrays
-- D_0, D_1, and D_2, and advanced recursively until the final
-- D_0 = SUM {C_k * Q_k} is found. This will be named Poly_Coefficients.
--
--     D_n+2(X) = 0
--     D_n+1(X) = 0
--     D_m(X)   = C_m + (X - Alpha_m+1)*D_m+1(X) - Beta_m+2*D_m+2(X)
--
-- where n = Desired_Poly_Degree and m is in 0..n.
-- The arrays D_0, D_1, and D_2 are already initialized.
-- They contain the coeff's of powers of X in the orthogonal
-- polynomials D_m, D_m+1, D_m+2.  The formulas above
-- imply:  for m in Desired_Poly_Degree .. 0:
--
--   D_0(k) := D_1(k-1) - Alpha(m)*D_1(k) - Beta(m)*D_2(k)
--   D_0(0) := D_0(0) + C(m)
--
-- where k goes from 1 to Desired_Poly_Degree in the first equation.
-- Notice that the above formula must be modified at m = n and m = n-1.
-- This gives the Poly coefficients for X in the range [-2,2].  Some
-- poly shifting has to be done to get the coefficients in the the
-- true range.  (Common source of floating point overflow.)
--
-- The algorithm can be generalized to functions F of many variables (2).
--
-- (1) Forsythe, G. E., J. SIAM 5 (1957), 74-87.
--     This package follows Shampine:
--     Shampine and Allen, Numerical Computing, An Introduction,(1973).
-- (2) Bartels and Jezioranski, ACM Transactions on Math. Software,
--     Vol. 11 (1985), p. 201-217.
--

generic

   type Real is private;
   --  The package can be instantiated with either conventional Floats,
   --  or with an abstract (extended precision) floating point.

   type Points_Index is range <>;
   type Data is array (Points_Index) of Real;

   Zero : Real;
   One  : Real;
   --  For initializing objects of private type Real. 'Real' can be
   --  an extended precision abstract float type.

   with function "*" (X : Real; Y : Real) return Real is <>;
   with function "+" (X : Real; Y : Real) return Real is <>;
   with function "-" (X : Real; Y : Real) return Real is <>;
   with function "/" (X : Real; Y : Real) return Real is <>;
   with function "<" (X : Real; Y : Real) return Boolean is <>;
   with function "**" (X : Real; N : Integer) return Real is <>;
   with function "-" (X : Real) return Real is <>;

package Orthogonal_Polys is

   type Integer32 is range -2**31 .. 2**31-1;

   Max_Order_Of_Poly : constant Integer32 := Integer32 (Data'Length - 1);
   --  May want to set this to a smaller value if the number of
   --  data points is >> the max order of polys.  In any case,
   --  this constant cannot be greater than Data'Length - 1.


   --  Section 1.
   --
   --  Data structures for polynomials.
   --  Type Poly_Data holds the actual values of individual Polys at X-axis
   --  grid points X_j.  Also holds a normalization factor for the Poly.
   --  Type Polynomial holds coefficients Alpha and Beta for generating
   --  Polys recursively via the Gram-Schmidt method.  (Plus other things.)

   subtype Coeff_Index is Integer32 range 0..Max_Order_Of_Poly;
   --  Index for array of polynomial coefficients.  The Max_Order must be
   --  less than the number of data points.  An assertion verifies this.
   --  (The max number if data points is taken to be Data'Length.)

   type Poly_Data is record
      Points  : Data         := (others => Zero);
      First   : Points_Index := Points_Index'First;
      Last    : Points_Index := Points_Index'Last;
      Squared : Real         := Zero;
      Degree  : Coeff_Index  := Coeff_Index'First;
   end record;
   --  Value of the polynomial at the X axis grid points.  Squared will
   --  contain a factor for normalizing the polys.  The Polys are not
   --  contained in normalized form in the field Points.  Squared will
   --  contain Inner_Product (Poly, Poly, First, Last, Weights).
   --  Poly_data.Points holds the explicit Y values of the poly defined at
   --  points X_j .. X_k, where j and k are Data_First amd Data_Last as
   --  stored in Points.  These end points are input by the user on
   --  a call to Start_Gram_Schmidt_Poly_Recursion.  This is where the
   --  Data_First and Data_Last, (start and end of the Polys) are input,
   --  nowhere else.  After this, all operations on this data are in
   --  the range Data_First..Data_Last.  To turn off data points in that
   --  range you set their associated weights to 0.0.

   type Recursion_Coeffs is array (Coeff_Index) of Real;
   --  The Alpha and Beta coefficients for generating polynomial recursively.
   --  These are calculated in calls to Get_Next_Poly and
   --  Start_Gram_Schmidt_Poly_Recursion.

   type Polynomials is limited private;
   --  Contains the Gram_Schmidt recursion coefficients Alpha and Beta
   --  for generating polynomials up to some order.  That order is given
   --  by Degree_Of_Poly, also contained in Polynomials.  The only way
   --  to use type Polynomials is to initialize it with a call to
   --  Start_Gram_Schmidt_Poly_Recursion, and to update it with calls
   --  to Get_Next_Poly.


   --  Section 2.
   --
   --  Construct the unique set of orthogonal polys associated with
   --  X-axis grid points "X_axis" and weights "Weights", on range
   --  First .. Last of Points_Index.

   procedure Start_Gram_Schmidt_Poly_Recursion
     (X_axis         : in     Data;
      Weights        : in     Data;
      First, Last    : in     Points_Index;
      Poly_0, Poly_1 : in out Poly_Data;
      Poly_Set       :    out Polynomials);

   --  Starts off the process of generating polys.  All computation
   --  is henceforth done on interval First..Last of Points_Index.
   --  If you want to leave out point (X(i), Y(i)) during the fitting,
   --  then set the corresponding weight to
   --          Weight(i) := Zero;
   --  Typical setting for Weights is simply Weights := (others => 1.0).
   --  This is the only place that First, Last, and X_axis are input. They
   --  are stored in Poly_Set and Poly_N, and subsequent operations
   --  read them from these objects.
   --  The in out parameters (Poly_O and Poly_1) need no initialization.

   procedure Get_Next_Poly
     (Poly_0, Poly_1 : in     Poly_Data;
      Weights        : in     Data;
      Poly_2         : in out Poly_Data;
      Poly_Set       : in out Polynomials);

   --  Generate successive polys given the 2 previous polys.
   --  This procedure creates Poly_2, and updates Poly_Set so that it
   --  includes recursion coefficients for generating Poly_2.
   --  "Weights" must contain the same weights used in the call to
   --  Start_Gram_Schmidt_Poly_Recursion that started off the
   --  process of creating poly_0 and poly_1.


   --  Section 3.
   --
   --  Section contains procedures that operate on type Polynomial.
   --  The procedures sum linear combinations of orthogonal polynomials, and
   --  calculate values, derivatives and integrals of these sums.
   --  The coefficients
   --  of the polynomials in the linear combinations, called C_k,
   --  will be stored in array C of type Poly_Sum_Coeffs.
   --  If you want the value of a single Poly, its derivatives, or
   --  integrals, then you use a delta-function for C_k: let C_k = 0.0
   --  everywhere except at the index of the desired Poly, which you set to
   --  1.0.  Can also get coefficients of powers of X, but be aware that
   --  this is an error prone process.


   function Max_Permissable_Degree_Of (P : Polynomials) return Coeff_Index;
   --  P is the unique set of (GS) orthogonal polys associated with a set
   --  of grid points {X_j} and a set of weights {W_j}.  The max degree
   --  of this set of orthogonal polys is limited by the number of distinct
   --  points in set {X_j}, (the number that have nonzero weight).  That
   --  number is returned by this function.

   subtype Poly_Sum_Coeffs is Recursion_Coeffs;
   --  Multiply these coefficients by polynomials, sum to make linear
   --  combinations.  These sums (usually Least squares
   --  fits) are evaluated by calls to Poly_Value, Derivatives_Of etc.

   function Poly_Value
     (X        : in Real;
      C        : in Poly_Sum_Coeffs;
      Poly_Set : in Polynomials)
      return Real;

   --  Returns the value of poly Q(X) = SUM {C_m * Q_m(X) }, at arbitrary X.
   --  The orthogonal polys Q_m are defined by Poly_Set.  Uses the Clenshaw
   --  summation formula.
   --  (The procedure Polynomial_Fit, by the way, only calculates
   --  Q for values of X that are in the data set (X_j), and returns them in
   --  the array Best_Fit_Poly.)

   function Poly_Integral
     (X                    : in  Real;
      C                    : in  Poly_Sum_Coeffs;
      Poly_Set             : in  Polynomials;
      Order_Of_Integration : in Coeff_Index := 1)
      return Real;

   --  Returns the indefinite integral of poly Q(X) = SUM [ C_k * Q_k(X) ].
   --  The orthogonal polys Q_k are defined by Poly_Set.  To get the
   --  integral of Q on [A, B] use Integral(B) - Integral(A).  Uses Clenshaw
   --  summation formula to get coefficients of powers of X, then integrates
   --  by summing the poly with Horner's rule.  Not the best way to do
   --  quadrature in general, but it's a good way to get integrals of
   --  polynomial least-square fits to data.

   subtype Derivative_List is Recursion_Coeffs;
   subtype Derivatives_Index is Coeff_Index; -- Starts at 0.

   procedure Poly_Derivatives
     (Derivatives    : in out Derivative_List;
      X              : in     Real;
      Order_Of_Deriv : in     Derivatives_Index;
      C              : in     Poly_Sum_Coeffs;
      Poly_Set       : in     Polynomials);

   --  Returns the derivative of the polynomial
   --  Q(X) = SUM [ C_k * Q_k(X) ].   The m-th derivative is returned
   --  as the m-th component of the array Derivatives.  All lower order
   --  derivatives are return also, since they must be calculated along
   --  along the way.  You may not need them, but when you do, it is
   --  better not to call the procedure m times (when you can call it
   --  only once).
   --
   --  Procedure Poly_Derivatives provides another way of
   --  of getting the coefficients (E) of powers of X.  The m-th derivative
   --  of the polynomial, evaluated at X = 0, equals E_m multiplied
   --  by m!.  In many ways this method may be superior to that used in
   --  in procedure Get_Coefficients_Of_Powers_Of_X.  It's likely
   --  to work in situations where Get_Coeffs... fails, (and vice-versa).
   --  Also, if you calculate a high degree polynomial fit, you may use
   --  Poly_Derivative to get just the first few E_m, rather than all
   --  of them, which is what Get_Coeffs_Of_Powers_Of_X does, and which
   --  may lead to floating point exceptions.

   type Powers_Of_X_Coeffs  is array (Coeff_Index) of Real;

   procedure Get_Coeffs_Of_Powers_Of_X
     (Coeffs   :    out Powers_Of_X_Coeffs;
      C        : in     Poly_Sum_Coeffs;
      Poly_Set : in     Polynomials);

   --  Often one wants the coefficients of the powers of X in the
   --  least-squares polynomial, instead of the coefficients C_m of the
   --  orthogonal polynomials Q_m.  For example, one commonly wants
   --  the best linear fit to the data, or the best parabolic fit
   --  to the data,  F = E_0 + E_1 * X + E_2 * X**2,
   --  and so one wants the coefficients E_m in the sum:
   --
   --                          N
   --   Best_Fit_Poly (X)  =  SUM ( E_m * X**m ).
   --                         m=0
   --
   --  The coefficients E_m are returned in the array Coeffs.
   --  This procedure attempts to get all of the coefficients of the powers
   --  of X up to the degree of the polynomial that is contained in
   --  the Data structure Poly_Set.  This polynomial degree is usually,
   --  but not always, the degree set by the parameter Desired_Poly_Degree
   --  in the procedure Poly_Fit.
   --  If you calculate a high order polynomial least squares fit,
   --  it may not be a good idea the use this routine. This is a common
   --  source of floating point overflow, catastrophic loss of precision,
   --  etc.  These problems depend on a lot of things, so you can give it
   --  a try if you want, but this routine is really mostly for
   --  getting the best fit line, or parabola.  In these cases set
   --  Desired_Poly_Degree to 1 or 2 in procedure polynomial_fit to get
   --  the polynomial fit Poly_Set.  Then send Poly_Set into this routine.
   --
   --  The procedure Poly_Derivatives above provides another way of
   --  of getting the coefficients of powers of X.  See the above note.

   function Inner_Product
     (X, Y        : in Data;
      First, Last : in Points_Index;
      Weights     : in Data)
      return Real;

   --  Usually you project the polynomials onto the data with the
   --  function Inner_Product in order to make least squares fits.
   --  Also useful in testing orthogonality of polys.  This is the
   --  weighted inner product over which the polynomials (calculated
   --  below) are orthogonal.

   procedure Poly_Fit
     (Data_To_Fit         : in     Data;
      X_axis, Weights     : in     Data;
      First, Last         : in     Points_Index;
      Desired_Poly_Degree : in     Coeff_Index;
      Best_Fit_Poly       : in out Poly_Data;
      Best_Fit_Coeffs     : in out Poly_Sum_Coeffs;
      Poly_Set            : in out Polynomials;
      Mean_Square_Error   :    out Real);

private

   Two  : constant Real := One + One;
   Four : constant Real := Two + Two;

   Smallest_Delta_X  : constant Real := Two ** (-128);
   --  Some arbitrary small number.  Another choice might be:
   --  Smallest_Delta_X  : constant Real := 2.0**(-Real'Safe_Emax / 2);
   --  If Data points are not separated by greater than Smallest_Delta_X,
   --  then an exception will be raised.
   --  Prevents some of the more immediate of the possible overflows.
   --  This number is comparable to the square root of Real'Safe_Small
   --  because Safe_Small is 2.0**(-Safe_Emax-1). (Safe_Emin = -Safe_Emax.)
   --  A typical Safe_Small for IEEE 754 is 10**(-308).  The number
   --  here would be approximately 10**(-154).

   Smallest_Weight : constant Real := Two ** (-256);
   --  Set it to some arbitrary small number.  Another choice
   --  might be Smallest_Weight : constant Real := Real'Safe_Small;
   --  If a weight is set below this value, then we pretend
   --  that the data point is not in the set.

   Max_No_Of_Data_Points : constant Integer32 := Integer32 (Data'Length);
   --  So can verify that polynomial degree is less than
   --  the number of data points.

   type X_Axis_Scale is record
      Slope : Real := One;
      Const : Real := Zero;
   end record;

   type Polynomials is record
      Alpha          : Recursion_Coeffs  := (others => Zero);
      Beta           : Recursion_Coeffs  := (others => Zero);
      Scale_Factors  : X_Axis_Scale;
      X_scaled       : Data              := (others => Zero);
      Degree_Of_Poly : Coeff_Index       := Coeff_Index'First;
      Max_Permissable_Degree_Of_Poly : Coeff_Index := Coeff_Index'First;
   end record;

   pragma Assert (Coeff_Index'First = 0);
   pragma Assert (Coeff_Index'Last >= 3);

end Orthogonal_Polys;

