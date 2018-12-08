
-----------------------------------------------------------------------
-- package body Crout_LU, LU decomposition, with equation solving
-- Copyright (C) 2008-2018 Jonathan S. Parker.
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

-- package Runge_Coeffs_pd_5
--
-- Package contains coefficients for the Cash-Karp
-- 5th order, 6 stage Runge Kutta method with error control.
--
-- A test routine is provided to make make sure that they
-- have not been corrupted.  If the coefficients do not pass 
-- the tests, then they are given in redundant and alternate 
-- forms at the end of this package (commented out).
--
-- NOTES ON THE ALGORITHM
--
-- Given a differential equ. dY/dt = F(t,Y) and an initial condition
-- Y(t), the Runge-Kutta formulas predict a new value of Y at
-- Y(t + h) with the formula
--
--                          max
--        Y(t + h) = Y(t) + SUM (K(j) * B(j))
--                          j=1
--
-- where max = Stages'Last, (which equals 13 here, since this is a 13
-- stage Runge-Kutta formula) where the quantities K are calculated from
--
--        K(1) = h * F(t, Y(t))
--
--                                         j-1
--        K(j) = h * F (t + C(j)*h, Y(t) + SUM (K(i)*A(j,i)) )
--                                         i=1
--
--  The Fehberg method, used here, provides two versions of the array
--  B, so that Y(t + h) may be calculated to both 7th order and to
--  8th order using the same K(i)'s.  This way, the local truncation
--  error can be estimated without calculating the K(i)'s twice.
--
--  If we apply this formula to a time-independent linear F, then we can
--  derive a condition that can be used to test the correctness of the
--  initialization values of the arrays A, B and C.  To derive such a
--  formula we use the fact that the RK prediction of Y(t + h) must equal
--  the Taylor series prediction up to the required order.  So,
--
--     K(1) = h*F*Y    (where F now is a matrix, and Y a vector)
--
--     K(2) = h*F*(Y + K(1)*A21)
--
--     K(3) = h*F*(Y + K(1)*A31 + K(2)*A32)
--
--     K(4) = h*F*(Y + K(1)*A41 + K(2)*A42 + K(3)*A43)
--
--  The linearity of F implies that F(a*Y + Z) = a*F*Y + F*Z so:
--
--     K(1) = h*F*Y
--
--     K(2) = h*F*Y +  A21*h^2*F^2*Y
--
--     K(3) = h*F*Y + (A31 + A32)*h^2*F^2*Y   +          A32*A21*h^3*F^3*Y
--
--     K(4) = h*F*Y + (A41+A42+A43)*h^2*F^2*Y + (A42*A21+A43*A31)h^3*F^3*Y
--                                              +       A43*A32*A21*h^4F^4*Y
--
--  Now we use the fact that we must have the RK prediction equal that
--  of the Taylor's series up to a certain order:
--
--     max
--     SUM (K(j) * B(j)) = h*F*Y + ... + (1/n!)*h^n*F^n*Y +  O(h^n+1)
--     j=1
--
-- Here n=8 for the coefficients B = B8 given below.  Its n=7 for B7.
-- The above formula gives us a relation between 1/n! and B and A.
-- This formula is used in the procedure TestRKPD given at the end
-- of this package.  We see immediately that we must have:
--
--  max
--  SUM (B(i))                  = 1/1!
--  i=1
--
--  max         i-1
--  SUM (B(i) * SUM(Aij)) = 1/2!
--  i=2         j=1
--

generic

   type Real is digits <>;

package Runge_Coeffs_pd_5 is

   subtype RK_Range  is Integer  range 0..6; 
   subtype Stages    is RK_Range range 0..6; -- always 0 .. 6
   type Coefficient  is array(RK_Range) of Real;
   type Coefficient_Array  is array(RK_Range) of Coefficient;


   procedure Test_Runge_Coeffs;


   A_rational : constant Coefficient_Array :=
  (
   (others => 0.0),
   (1.0/5.0,         others => 0.0),
   (3.0/40.0,        9.0/40.0,         others => 0.0),
   (44.0/45.0,      -56.0/15.0,       32.0/9.0,        others => 0.0),
   (19372.0/6561.0, -25360.0/2187.0,  64448.0/6561.0, -212.0/729.0, others => 0.0),
   (9017.0/3168.0,  -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0.0, 0.0),
   (35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0,   -2187.0/6784.0, 11.0/84.0, 0.0)
  );

   -- 4rth order:

   B4_rational : constant Coefficient :=
   (
    5179.0/57600.0,
    0.0,
    7571.0/16695.0,
    393.0/640.0,
   -92097.0/339200.0,
    187.0/2100.0,
    1.0/40.0
   );


   -- 5th order:

   B5_rational : constant Coefficient :=
   (
    35.0/384.0,
    0.0,
    500.0/1113.0,
    125.0/192.0,
   -2187.0/6784.0,
    11.0/84.0,
    0.0
   );

    
   --  coefficients C for getting Dt

   C_rational : constant Coefficient :=
   (
    0.0,
    1.0/5.0,
    3.0/10.0,
    4.0/5.0,
    8.0/9.0,
    1.0,
    1.0
   );

   C  : Coefficient renames  C_rational;
   B4 : Coefficient renames B4_rational;
   B5 : Coefficient renames B5_rational;
   A  : Coefficient_Array renames A_rational;

end Runge_Coeffs_pd_5;
