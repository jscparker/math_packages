
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
-----------------------------------------------------------------------

-- PACKAGE runge_pc_2
--
-- Exclusively for initializing Predictor-Corrector arrays.
--
-- The 8th order Runge Kutta formula of Prince and Dormand
-- (Journal of Computational and Applied Math.
-- vol. 7, p. 67 (1981)).
--
-- The program integrates     (d/dt)**2 Y = F (t, Y)    where t is the
-- independent variable (e.g. time), vector Y is the
-- dependent variable, and F is some function.

generic

   type Real is digits <>;
   --  The independent variable.  Its called Time
   --  throughout the package as though the differential
   --  equation dY/dt = F (t, Y) were time dependent.  Of cource
   --  it can have any meaning.

   type Dyn_Index is range <>;

   type Dynamical_Variable is array(Dyn_Index) of Real;
   --  The dependent variable.
   --  We require its exact form here so that some inner loop
   --  optimizations can be performed below.
   --  A two dimension array makes this package useful
   --  for a very wide range of problems.  If your Dynamical
   --  variable is really just a one-dimensional array, then you can still
   --  use it here (and with almost no loss in performance.)
   --  Just give the Dyn_Index a range 0..0.  ( Its best to do
   --  this to the Vector_Range_1 rather than Vector_Range_2.)
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
   --  operation must exist by the definition of the derivative,
   --  dY/dt = (Y(t+dt) - Y(t)) / dt

package runge_pc_2 is

  procedure Integrate
    (Final_Y              :    out Dynamical_Variable;
     Final_deriv_Of_Y     :    out Dynamical_Variable;
     Final_Time           : in     Real;
     Initial_Y            : in     Dynamical_Variable;
     Initial_deriv_Of_Y   : in     Dynamical_Variable;
     Initial_Time         : in     Real;
     No_Of_Steps          : in     Real);

end runge_pc_2;
