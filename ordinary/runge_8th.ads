
-----------------------------------------------------------------------
-- package Runge_8th, 8th order Prince and Dormand Runge-Kutta
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
-----------------------------------------------------------------------

-- PACKAGE Runge_8th
--
-- Package implements the 8th order Runge Kutta formula of Prince
-- and Dormand (Journal of Computational and Applied Math.
-- vol. 7, p. 67 (1981)).
--
-- Procedure Integrate solves  dY/dt = F (t, Y)  where t is the
-- independent variable (e.g. time), vector Y is the dependent 
-- variable, and F is some function. 
--
-- Just a reminder of why high order is important if you want high
-- accuracy: compare a 2nd order method (error per step ~ B*dt^3)
-- with a 9th order method (error per step ~ A*dt^10). Suppose each
-- method gives for this 'dt' a relative error of 10^(-6), (an error 
-- of 1 part per million) and we want to reduce this to 10^(-15).  
-- In other words we want to reduce the error by a factor of 10^(-9). 
-- If the 2nd order step size dt is divided by 2^10, then its error
-- (per step( falls by a factor (1/2^10)^3 = 1/2^30 = 1/10^9 as required.
-- If the 9th order step size dt is divided by 2^3, then its error 
-- falls by a factor (1/2^3)^10 = 1/2^30 = 1/10^9 as required. So 
-- the 2nd order method needs at least 128 times as many time steps
-- to perform an integration over the same time interval as the 9th
-- order method.  Even if the 2nd order method runs 5 times faster
-- per step than the 9th order method, then we should expect the
-- 9th order method to be ~25 times faster than the 2nd order method.
-- (Actually, globally the error term behaves as though it is an order
-- lower than the local (per step) error term, so the higher order is
-- even more favorable than the above argument indicates.)
--
-- Additional note: the error term per step of the Runge-Kutta looks
-- more like A*dt^10 than A*dt^9 as you near machine precision. (Prince
-- and Dormand optimized the coefficients to minimize this term.) Yet
-- another reason why this Runge-Kutta performs better than expected
-- in so many cases.
--
-- NOTES ON USE
--
-- Runge_8th assumes that the Dynamical Variable is a 1 dimensional 
-- array of Real numbers.  
--
-- To use this routine higher order differential equations have to
-- be reduced to 1st order equations.
-- For example a 3rd order equation in X becomes a first order
-- equation in the vector Y = (Y1,Y2,Y3) = (X, dX/dt, d/dt(dX/dt)).
-- If the equation in X is  (d/dt)**3 X = G(t,X), then the
-- equation in vector Y is
--
--           d/dt (Y1, Y2, Y3) = (Y2, Y3, G(t,Y1)).
--
-- If the equation in X is  d/dt(d/dt(dX/dt))) = G(t,X,dX/dt),
-- then the equation in Y is
--
--         d/dt (Y1, Y2, Y3) = (Y2, Y3, G(t,Y1,Y2)).
--
-- So the routine solves the equation dY/dt = F(t,Y), or, more
-- explicitly, d/dt (Y1, Y2, Y3) = (Y2, Y3, G(t,Y1,Y2)) = F(t,Y).
-- The user plugs in the function F(t,Y) as a generic formal function.
-- Even if F is t or Y independent, it must be in the form F(t,Y).
-- The user has to do all of the work to convert higher order
-- equations to first order equations. On some compilers,
-- performance is very sensitive to how this is done. That's
-- why this part is best left to the user.  I've seen factors of
-- 3 improvement in speed when this part of the data structure
-- is optimized.  A pragma Inline is also worth trying.
-- The generic formal parameters "+" and "*" below should be
-- compiled Inline.
--
-- Uses the elementary math functions Log and Exp (used in
-- the error control algorithm to implement the exponentiation 
-- function (**)). They're used each step, which can be slow!
--
-- Runge-Kutta routines can be used to integrate moderately
-- stiff equations: use constant step size without error control
-- and make the step size sufficiently small. (Routines specially
-- designed for stiff equations may be orders of magnitude faster.)
--
-- There's no fool proof way of getting the right answer when
-- you numerically integrate. In some cases error control does
-- not work well. In other cases, integration without
-- error control fails simply because it is too ineffient.
-- Using error control: if you have 15 digit floating point, 
-- then don't demand more than 12 or 13 digits of accuracy.
-- The error control algorithm begins to fail in that limit.

generic

   type Real is digits <>;
   --  The independent variable has type Real.  Its called Time
   --  throughout the package as though the differential
   --  equation dY/dt = F (t, Y) were time dependent.

   type Dyn_Index is range <>;

   type Dynamical_Variable is array(Dyn_Index) of Real;
   --  The dependent variable.  Declared as an array here so that
   --  some inner loop optimizations can be performed below.
   --
   --  To use this routine, reduce higher order differential
   --  equations to 1st order.
   --  For example a 3rd order equation in X becomes a first order
   --  equation in the vector Y = (X, dX/dt, d/dt(dX/dt)).

   with function F
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines the equation to be integrated: dY/dt = F (t, Y).  Even if
   --  the equation is t or Y independent, it must be entered in this form.

   with function "*"
     (Left  : Real;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Multiplication of the independent by the dependent variable. 
   --  An operation of this sort exists by definition of the
   --  derivative,  dY/dt = (Y(t+dt) - Y(t)) / dt,
   --  which requires multiplication of the (inverse) independent
   --  variable t, by the Dynamical variable Y (the dependent variable).

   with function "+"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Summation of the dependent variable.
   --  The operation must exist by the definition of the derivative,
   --  dY/dt = (Y(t+dt) - Y(t)) / dt =  (Y(t+dt) + (-1)*Y(t)) / dt.

   with function "-"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;

   with function Norm
     (Y1: in Dynamical_Variable) return Real is <>;
   --  For error control we need to know the distance between two objects.
   --  Norm define a distance between the two vector, Y1 and Y2
   --  using: distance = Norm (Y1 - Y2),
   --  For example, if Dynamical_Variable is complex, then Norm can 
   --  be  Sqrt (Real (Y1(1)) * Conjugate (Y1(1))));
   --  If it is a vector then the metric may be
   --  Sqrt (Transpose(Y1(1)) * Y1(1)), etc.
   --
   --  Recommended:  Sum the absolute values of the components of the vector.

package Runge_8th is

   type Step_Integer is range 1 .. 2**31-1;

   --  The step size Delta_t is never input. Instead
   --  you input Starting_Time, Final_Time, and No_Of_Steps.
   --
   --  If error control is disabled, then Delta_t is
   --  (Final_Time - Starting_Time) / No_Of_Steps. 
   --
   --  If error control is enabled, then step size Delta_t is variable,
   --  and procedure Integrate chooses its own Delta_t each time step.
   --  Each time step, Delta_t is chosen to reduce fractional error
   --  per step to something less than Error_Tolerance. It does not
   --  and cannot reduce global error to any given value.
   --
   --  All components of vector Initial_State must be initialized.


   procedure Integrate
     (Final_State           :    out Dynamical_Variable;
      Final_Time            : in     Real;
      Initial_State         : in     Dynamical_Variable;
      Initial_Time          : in     Real;
      No_Of_Steps           : in     Step_Integer;
      Error_Control_Desired : in     Boolean       := False;
      Error_Tolerance       : in     Real          := +1.0E-10);

private

   Zero : constant := +0.0;
   Half : constant := +0.5;
   One  : constant := +1.0;
   Two  : constant := +2.0;
   Four : constant := +4.0;
   Nine_Tenths   : constant := +0.9;
   One_Eighth    : constant := +0.125;

   Test_of_Runge_Coeffs_Desired : constant Boolean := False;
   --  The test (exported by Runge_Coeffs_13) prints an error message if
   --  failure is detected. Should be set to True during testing.

end Runge_8th;
