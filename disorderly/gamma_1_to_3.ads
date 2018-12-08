
--
--  The following is based on John Maddock's extended precision
--  (Boost Library) gamma function.  It uses the gamma rational
--  poly functions (for the range x = 1 to 3) worked out by John Maddock,
--  so his Copyright will be retained.
--
-------------------------------------------------------------------------- 
--  (C) Copyright John Maddock 2006.
--  Use, modification and distribution are subject to the
--  Boost Software License, Version 1.0. (See accompanying file
--  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
--
--Boost Software License - Version 1.0 - August 17th, 2003
--
--Permission is hereby granted, free of charge, to any person or organization
--obtaining a copy of the software and accompanying documentation covered by
--this license (the "Software") to use, reproduce, display, distribute,
--execute, and transmit the Software, and to prepare derivative works of the
--Software, and to permit third-parties to whom the Software is furnished to
--do so, all subject to the following:
--
--The copyright notices in the Software and this entire statement, including
--the above license grant, this restriction and the following disclaimer,
--must be included in all copies of the Software, in whole or in part, and
--all derivative works of the Software, unless such copies or derivative
--works are solely in the form of machine-executable object code generated by
--a source language processor.
--
--THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
--IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
--FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
--SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
--FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
--ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
--DEALINGS IN THE SOFTWARE.
-------------------------------------------------------------------------- 
-- 
-- Natural logarithm of Gamma function for positive real arguments.
--
generic

   type Real is digits <>;

package Gamma_1_to_3 is

   pragma Pure(Gamma_1_to_3);

   function Log_Gamma_1_to_3 (x : in Real) return Real;
   -- Only good for     1 < x < 3

private

   Real_Epsilon : constant Real := (+0.125) * Real'Epsilon;
   -- have to modify this if Real is abstract

end Gamma_1_to_3;

