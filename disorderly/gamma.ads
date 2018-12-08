
-------------------------------------------------------------------------------
-- package Gamma, Log of Gamma Function
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

-- package Gamma
-- 
-- Natural logarithm of Gamma function for positive real arguments.
--
-- Uses Stieltjes' Continued Fraction method arg > 14, and rational
-- polynomial approximations for 0 < arg < 16.
--
-- Meant to be good to better than 30 digits (if you can instantiate
-- with a 30 digit Real). That's why the rational
-- polynomial approximations are so complicated at x < 16, and
-- why the order of the Stieltjes' Continued Fraction is so high.
-- But ordinarily you use a 15 digit Real.

generic

   type Real is digits <>;

package Gamma is

   pragma Pure (Gamma);

   function Log_Gamma (x : in Real) return Real;
   -- For x > 0 only.  Natural logarithm of Gamma function.
   --
   -- Uses Stieltjes' Continued Fraction method above x=14.
   -- For 0 < x < 10 uses simplest rational approx.
   -- Probably good to a lot better than 32 digits. 
   
   function Log_Gamma_0_to_16 (x : in Real) return Real;
   -- This is the part that uses a rational polynomial approx.
   -- Only good for   0 < x < 16   !

   procedure Test_Stieltjes_Coefficients;

private

   Real_Epsilon : constant Real := (+0.125) * Real'Epsilon;
   -- have to modify this if Real is abstract extended

end Gamma;

