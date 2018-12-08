
-- PACKAGE e_Derivs
--
-- Automatic differentiation routines for calculating arbitrary order
-- derivatives of arithmetical expressions constructed out of operators
-- "+", "-", "*", "/", "**", Compose, Sin, Cos, Exp, Log, Reciprocal, 
-- Sqrt, Arccos, Arsin, Arctan, or any user defined function.  Also
-- includes routines for summing the Taylor series for the value of the
-- Function,  its Integrals, and Derivatives at arguments near to those
-- at which the original Function, and its derivs were evaluated.
-- 
-- Warnings: 
-- High order differentiation is prone to overflows/underflows and
-- loss of accuracy, especially near singularities. For example 
-- the high order derivs of 1/t, log(t) or Sqrt(t) go as high powers
-- of 1/t, so if t is small or large enough, then exceptions happen.
-- Also note that numerical accuracy can degrade rapidly with order of 
-- derivative.  Nothing clever has been done to improve numerical accuracy.
-- If you use 15 digit floats, (with Real'Epsilon around 1.0e-15), then error
-- many times this magnitude accumulates, so that multiplying by, say,
-- 16! to get Un-reduced (true) derivatives can destroy all accuracy.
--
-- Results are returned as Leibniz Reduced Derivatives (n-th derivative
-- divided by n!), which are usually better behaved numerically than
-- the high order Un-reduced derivatives.
--
-- Here is how it works: you place the first N derivatives of some
-- function F(t) in array F, and the first N derivatives of another function
-- G(t) in array G, then the operator "*" will calculate the 1st N 
-- derivatives of the product F(t)*G(t) and place the result into array H when
-- you write H := F * G.  Similarly, H := Compose (F, G) places the first N
-- derivatives of the composition of the 2 functions F(G(t)) into array H,
-- provided F is the set of derivatives of F evaluated at G(t).  Routines
-- are provided to help calculate high order derivatives of elementary
-- functions such as Sin and Log.  (The routines have names like Sin_d and
-- Log_d, respectively.)
-- You do not have to use the functions provided.  If you determine the
-- derivatives of a function numerically, (for example a function that is the
-- solution of a numerically integrated differential equ.) and place them in an
-- array of type Derivatives, then this works just as well.
--
-- The program does not parse expressions for you.  You must translate
-- expressions like F(G(t)) into: Compose (F, G).  However, most other steps in
-- the translations use the original expression in unmodified form. Use
-- Parentheses to indicate precedence of the operators.
--
-- The derivatives are input and output in Leibniz Reduced form.
-- (To get Leibnitz reduced derivatives you divide the N-th derivative by N!.)
-- All computation is done internally in reduced form,
-- because this is usually the most efficient way.  You can translate 
-- between Reduced and Unreduced forms of derivatives with procedures
-- Make_Reduced, and Un_Reduce.
-- 
-- Notice that the reduced derivatives of H(t) are also the Taylor coefficients
-- of function H(t).  More precisely, the Taylor's series for H(t), (if we know
-- (d/dt)**N H(t_0), the set of derivatives of H evaluated at point t_0), is
-- the sum over coefficents (t - t_0)**N * (d/dt)**N H(t_0) / N!. The functions
-- in this package can be used to directly calculate (d/dt)**N H(t_0) / N!,
-- if H is composed of functions whose derivatives are known.  For example,
-- if we know the derivatives of F and of G, then we can calculate Taylor
-- series for H(t) = F(G(t)).  So the routines in this package might be thought
-- of as routines for gettin Taylor series coefficients of complicated functions, 
-- given Taylor Series of the simpler functions from which they are constructed.
-- For example suppose you numerically integrate 2 differential equations for
-- F and G, returning the high order derivatives of F and G along the way.
-- Then if you want (say) the integral of F*G, you calculate the Taylor
-- series of H = F*G at various points t_0, and then calculate the area under
-- the Taylor Series polynomial by elementary means.  Similarly, you can
-- numerically integrate differential equations by the Taylor series method
-- if these routines are used to calculate the derivatives of the
-- function that defines the differential equation.
--
-- Exceptions:
-- Some functions raise d_Argument_Error if arguments are out of range.
--
generic 

   Max_Order_Of_Deriv : Positive;

   type Real is private;

   type Real_8 is digits <>;

   --  Never have to enter these explicitly:  
   with function Sin (X : Real) return Real is <>;
   with function Cos (X : Real) return Real is <>;
   with function Exp (X : Real) return Real is <>;
   with function Log (X : Real) return Real is <>; -- base e
   with function Sqrt (X : Real) return Real is <>; 
   with function Arcsin (X : Real) return Real is <>; 
   with function Arccos (X : Real) return Real is <>; 
   with function Arctan (X : Real) return Real is <>; 

   with function "+" (X : Real_8) return Real is <>; 
   -- use for translating 15 digit float  (Real_8) to extended (Real).
   -- used here only for making One_d, Half, etc.

   with function "-" (X : Real) return Real is <>;
   with function "+" (X, Y : Real) return Real is <>;
   with function "-" (X, Y : Real) return Real is <>;
   with function "*" (X, Y : Real) return Real is <>;
   with function "/" (X, Y : Real) return Real is <>;
   with function "<=" (X, Y : Real) return Boolean is <>;
   with function ">=" (X, Y : Real) return Boolean is <>;
   with function  "=" (X, Y : Real) return Boolean is <>;
   with function "**" (X : Real; Exponent : Natural) return Real is <>;
   --need Equality test for Reciprocal.
  
package e_Derivs is
 
   One_d  : constant Real := (+1.0);
   Zero_d : constant Real := (+0.0);

   subtype Deriv_Index is Integer range 0..Max_Order_Of_Deriv;
   
   type Derivatives is array(Deriv_Index) of Real;

   function "+" (F, G : in  Derivatives) return Derivatives;
   --  Returns reduced derivatives of F(t) + G(t), given reduced derivs of F, G.

   function "-" (F, G : in  Derivatives) return Derivatives;
   --  Returns reduced derivatives of F(t) - G(t), given reduced derivs of F, G.

   function "*" (F, G : in  Derivatives) return Derivatives;
   --  Returns reduced derivatives of F(t) * G(t), given reduced derivs of F, G.

   function "/" (F, G : in  Derivatives) return Derivatives;
   --  Returns reduced derivatives of F(t) / G(t), given reduced derivs of F, G.

   function "**" 
     (F        : in Derivatives;
      Exponent : in Natural) 
      return Derivatives;
        
   function Compose (F, G : in Derivatives) return Derivatives;
   --  Returns reduced derivatives of F(G(t)), given reduced derivs of F, G.
   --  F must hold the reduced derivs of F(t'), evaluated at t'=G(t).  G holds
   --  reduced derivs of G(t) pre-evaluated at whatever t the user desires.
   --  So the derivatives of F(G(t)) are returned by: Compose(F(t'),G).
   
   function Reciprocal (F : in  Derivatives) return Derivatives;
   --  Returns reduced derivatives of 1 / F(t), given reduced derivs of F.
     
   procedure Make_Reduced (F : in out Derivatives);
   --  Make reduced derivatives out of ordinary derivs, by dividing by Order!.

   procedure Un_Reduce (F : in out Derivatives);
   --  Make ordinary derivatives out of reduced derivs, by multiplying by Order!
   
   d_Argument_Error : exception;
   --  A little checking of arguments is done, but most is left to the
   --  functions of type Real that the procedures here are built out of.

   
   --**************************************************************************
   -- Use taylor series to evaluate the Function, its Integrals and Derivatives
   -- at argument t_1 /= t_0 (where t_0 is the pt. at which the original Function
   -- and its derivs were evaluated).
   --**************************************************************************

   function Taylor_Sum
     (t                : in Real;
      F                : in Derivatives;
      Max_Index        : in Deriv_Index   := Deriv_Index'Last;
      No_Of_Iterations : in Natural       := 0)
      return Real;
   --  Given F and reduced derivs of F at t_0, routine sums taylor series for
   --  F(t_1), where t_1 = t + t_0. (Summing polys is always highly error prone.)

   function Taylor_Sum_2
     (t                : in Real;
      F                : in Derivatives;
      Max_Index        : in Deriv_Index   := Deriv_Index'Last) return Real;
   --  Slightly different  algor.

   function Derivative_Of 
     (F              : in  Derivatives;
      Order_Of_Deriv : in Deriv_Index)
      return Derivatives;
   --  Returns reduced derivatives of (d/dt)^k F(t), given reduced derivs of F,
   --  where k = Order_Of_Deriv.
     
   function Integral_Of 
     (F : in  Derivatives)
      return Derivatives;
   --  Integration const is set to 0.  To get area under curve F(t),
   --  say in [t_0 - dt, t_0 + dt], let F be reduced derivs of F(t) at t_0,
   --  area = Taylor_Sum(dt, Integral_Of(F)) - Taylor_Sum(-dt, Integral_Of(F)).
     
   --**************************************************************************
   -- The following functions operate on argument of type Real rather than type
   -- Derivative.  They return the high order *reduced* derivatives of common
   -- functions to make it easier to construct more complicated functions.
   -- Their arguments are first order polynomials of "time", also for convenience.
   -- I.e., Sin_d returns reduced derivs of Sin (constant1 * time + constant2).
   -- Could use the Composition operator, but that would be inefficient.
   --**************************************************************************

   function "**"
     (Time     : in Real;
      Exponent : in Natural) return Derivatives;
        
   function Reciprocal (Time : in  Real) return Derivatives;

   function Sin_d 
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d) return Derivatives;
       
   function Cos_d 
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d) return Derivatives;
       
   function Exp_d
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d) return Derivatives;
   
   function Log_d
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d) return Derivatives;
   
   function Sqrt_d
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d) return Derivatives;
   
   function Arctan_d
     (Time : in Real) return Derivatives;
   
   function Arcsin_d
     (Time : in Real) return Derivatives;
   
   function Arccos_d
     (Time : in Real) return Derivatives;  -- TEST THIS.
   
end e_Derivs;
