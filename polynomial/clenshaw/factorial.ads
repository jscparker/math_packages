
-- PACKAGE Factorial
-- 
-- Natural logarithm of Factorial for arguments 0 and greater.
--
-- Uses Stieltjes' Continued Fraction method arg > 32, and
-- table for 0 <= arg <= 32.
--
-- Meant to be good to (much much) better than 30 digits (if you
-- can instantiate with a 30+ digit Real). That's why the
-- order of the Stieltjes' Continued Fraction is so high.
-- But ordinarily you use a 15 or 18 digit Real.
--
generic

   type Real is digits <>;

package Factorial is

   pragma Pure (Factorial);

   function Log_Factorial (N : in Natural) return Real;
   -- For N >= 0.  Natural logarithm of N!
   --
   -- Uses Stieltjes' Continued Fraction method above N=20.
   -- Probably good to a lot better than 32 digits (if you can find
   -- that many digits).
   
   procedure Test_Stieltjes_Coefficients;
   procedure Test_Log_Factorial_Table;
   --  Should call these sometime. The test routine will do it for you.
   --  They raise Program_Error if the tables of constants have mutated.

end Factorial;

