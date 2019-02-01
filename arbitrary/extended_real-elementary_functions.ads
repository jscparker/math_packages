
-- PACKAGE Extended_Real.Elementary_Functions
--
-- Taylor series are used to get Exp, Sin and Cos.  Once these are calculated,
-- Newton-Raphson iteration gives inverse functions: Log, Arccos, and Arcsin,
-- Arctan.  Similarly, starting with the function G(Y) = Y**(-N),
-- Newton-Raphson iteration gives the inverse: F(X) = X**(-1/N),
-- the reciprocal of the Nth root of X.  Newton-Raphson is used directly to get
-- Sqrt, Inverse, and Inverse_Nth_root. Inverse(X) is the reciprocal of
-- of X.  It's usually faster than the "/" function, (One / X).
--
-- Newton-Raphson iteration is used to calculate Y = F(a) when F's inverse
-- function G(Y) is known.  (G satisfies G(F(a)) = a.)  Say we want Y = F(a),
-- We can't calculate F, but given a Y we can calculate G(Y) and we know a. Then
-- use Newton-Raphson to solve for Y in equation G(Y) = a.  The iteration is
-- 
--      dG/dY(Y_0) = (G(Y_0) - a) / (Y_0 - Y_1), or,
--
--      Y_1 = Y_0 - (G(Y_0) - a) / dG/dY(Y_0).
--
-- For example, if want F(a) = Log(a), and we can get G(Y) = Exp(Y), then we
-- iterate for Y = Log(a) using:
--
--      Y_1 = Y_0 - (Exp(Y_0) - a) / Exp(Y_0).
--
-- Similarly, if we want F(a) = a**(-1/N) and we know G(Y) = Y**(-N), then:
--
--      Y_1 = Y_0 - (Y_0**(-N) - a) / (-N*Y**(-N-1),
--
--          = Y_0 - (1 - a*Y_0**N) * Y / (-N).
--
-- Argument reduction is necessary in most of the routines.  Some of the
-- arg reduction ideas come from Brent, D M Smith, D H Bailey. 
-- 

generic

   --  Functions for type Real.  Must be correct to 48 bits.
   --  You should never have to enter these parameters, as long as
   --  they are visible at instantiation.   Log is natural.
   
   with function   Sqrt (X : Real) return Real is <>;
   with function    Log (X : Real) return Real is <>; -- Natural log. (ln)
   with function    Exp (X : Real) return Real is <>; -- inverse of log
   with function Arcsin (X : Real) return Real is <>; -- inverse of sin
   with function Arctan (Y : Real; X : Real := 1.0) return Real is <>; 
   
package Extended_Real.Elementary_Functions is

   function   Sqrt (X : e_Real) return e_Real;
   function    Exp (X : e_Real) return e_Real;
   function    Log (X : e_Real) return e_Real;  
   function    Log (X : e_Real; Base : e_Real) return e_Real;  
   function    Sin (X : e_Real) return e_Real;
   function    Sin (X : e_Real; Cycle : e_Real) return e_Real;
   function    Cos (X : e_Real) return e_Real;
   function    Cos (X : e_Real; Cycle : e_Real) return e_Real;
   function Arcsin (X : e_Real) return e_Real;
   function Arccos (X : e_Real) return e_Real;
   function   "**" (Left : e_Real; Right : e_Real) return e_Real;  

   function Arctan (X : e_Real) return e_Real;
   --  Output is in range [-Pi/2 .. Pi/2] only. Arctan (Infinity) = Pi/2.

   function Reciprocal (X : e_Real) return e_Real;
   --  Newton-Raphson inversion.  Usually faster than One / X.

   function Divide (Z, X : e_Real) return e_Real;
   --  Newton-Raphson Z / X.

     
   function Reciprocal_Nth_Root (X : e_Real; N : Positive) return e_Real;
   --  Reciprocal of the N-th root of X:  X**(-1/N) = 1 / X**(1/N).  One way
   --  to get the N-th root of X is to take One / Reciprocal_Nth_Root(X, N).
   --  N must be less than Radix - 1, which is usually 2**29 - 1.
   --  (This function is non-standard, but is used by some of the other
   --   routines, so might as well export it.)
   
   function Reciprocal_Sqrt (X : e_Real) return e_Real;

   
   function e_Quarter_Pi         return e_Real;  -- returns Pi/4 by arcsin method.
   function e_Log_2              return e_Real;  -- returns Log(2.0).
   function e_Inverse_Pi         return e_Real;  -- returns 1/Pi by arcsin method.
   function e_Inverse_Sqrt_2     return e_Real;  -- returns 1/Sqrt(2.0).
   function e_Half_Inverse_Log_2 return e_Real;  -- returns 0.5/Log(2.0).
   --  The above constants are calculated to max available precision.
   --  They are all slightly less than one - the the highest precision
   --  that this package is capable of.  So the above versions are
   --  preferred to the ones given below, which are somewhat
   --  greater than one.
   
   function e_Pi return e_Real;  -- returns Pi by arcsin method.
   
   E_Argument_Error : Exception;
   
end Extended_Real.Elementary_Functions;
