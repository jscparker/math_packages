
-----------------------------------------------------------------------
-- package Predictor_1, 17th order Predictor-Corrector integrator
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
------------------------------------------------------------------------

-- PACKAGE Predictor_1
--
-- Procedure Integrate is a 17th order Predictor-Corrector integration
-- routine for solving  dY/dt = F (t, Y). The independent variable is
-- t (e.g. time), vector Y is the dependent variable, and F is 
-- some function.
--
-- The Predictor-Corrector method uses least-squares to fit a 17th 
-- order polynomial to the 33 previous values of F in dY/dt = F(t,Y).
-- Least-squares-fit polynomials are then used to predict and 
-- correct the next values of F and Y.
--
-- NOTES ON USE
--
-- The following version assumes that the Dynamical Variable is
-- a 1 dimensional array of Real numbers.
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

generic

   type Real is digits <>;
   --  The independent variable.  It's called Time
   --  throughout the package as though the differential
   --  equation dY/dt = F (t, Y) were time dependent.  Of cource
   --  it can have any meaning.

   type Dyn_Index is range <>;

   type Dynamical_Variable is array(Dyn_Index) of Real;
   --  The dependent variable.
   --  We require its exact form here so that some inner loop 
   --  optimizations can be performed below.
   --
   --  To use this routine one must reduce higher order differential
   --  Equations to 1st order.
   --  For example a 3rd order equation in X becomes a first order
   --  equation in the vector Y = (X, dX/dt, d/dt(dX/dt)).

   with function F
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines the equation to be integrated as
   --  dY/dt = F (t, Y).  Even if the equation is t or Y
   --  independent, it must be entered in this form.

   with function "*"
     (Left  : Real;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines multiplication of the independent by the
   --  dependent variable.       An operation of this sort must exist by
   --  definition of the derivative, dY/dt = (Y(t+dt) - Y(t)) / dt,
   --  which requires multiplication of the (inverse) independent
   --  variable t, by the Dynamical variable Y (the dependent variable).

   with function "+"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines summation of the dependent variable.  Again, this
   --  operation must exist by the definition of the derivative,
   --  dY/dt = (Y(t+dt) - Y(t)) / dt =  (Y(t+dt) + (-1)*Y(t)) / dt.

   with function "-"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;

   with function Norm
     (Y1: in Dynamical_Variable) return Real is <>; -- not used yet.

package Predictor_1 is

   type Step_Integer is range 1 .. 2**31-1;

   procedure Integrate
     (Final_State   :    out Dynamical_Variable;
      Final_Time    : in     Real;
      Initial_State : in     Dynamical_Variable;
      Initial_Time  : in     Real;
      No_Of_Steps   : in     Step_Integer);
  
   --  You must init. all elements of array Initial_State.
   --
   --  The 1st 32 steps use an 8th order Runge-Kutta to initialize
   --  an array, so if No_Of_Steps < 33 then you are using the
   --  Runge-Kutta to do the integration.

  
   No_Of_Corrector_Evaluate_Iterations : constant Positive := 1;
   --  Notice there has to be at least 1 Correction step. That means
   --  this is always a predictor-corrector method. 1 usually best here.
  
   Final_Corrector_Desired : constant Boolean := True;
  
   Extrapolation_Desired : constant Boolean := True;

end Predictor_1;
