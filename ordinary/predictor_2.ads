--
-----------------------------------------------------------------------
-- package Predictor_2, 20th order Predictor-Corrector
-- Copyright (C) 2008-2009 Jonathan S. Parker.
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

-- PACKAGE Predictor_2
--
-- 20th order Predictor-Corrector for conservative differential equations.
--
-- The program integrates     (d/dt)**2 Y = F (t, Y)    where t is the
-- independent variable (e.g. time), vector Y is the
-- dependent variable, and F is some function. Certain higher order
-- equations can be reduced to this form.
--
-- Notice that F in independent of dY/dt.
--
-- NOTES ON USE
--
-- The following version assumes that the Dynamical Variable is
-- a 1 dimensional array of Real numbers.
--
-- The user plugs in the function F(t,Y) as a generic formal function.
-- Even if F is t or Y independent, it must be in the form F(t,Y).
-- The user must do all of the work to convert higher order
-- equations to 2nd order equations. 
-- 

generic

   type Real is digits <>;
   --  The independent variable.  Its called Time
   --  throughout the package as though the differential
   --  equation dY/dt = F (t, Y) were time dependent.

   type Dyn_Index is range <>;

   type Dynamical_Variable is array(Dyn_Index) of Real;
   --  The dependent variable.
   --  We require its exact form here so that some inner loop 
   --  optimizations can be performed below.

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
   --  Operation must exist by the definition of the derivative,
   --  dY/dt = (Y(t+dt) - Y(t)) / dt 

package Predictor_2 is

   --  You must init. all elements of array Initial_Y, Initial_deriv_Of_Y.

   procedure Integrate
     (Final_Y              :    out Dynamical_Variable;
      Final_deriv_Of_Y     :    out Dynamical_Variable;
      Final_Time           : in     Real;
      Initial_Y            : in     Dynamical_Variable;
      Initial_deriv_Of_Y   : in     Dynamical_Variable;
      Initial_Time         : in     Real;
      No_Of_Steps          : in     Real);


   No_Of_Corrector_Evaluate_Iterations : constant Integer := 1;
   --  1 is good here.

   Final_Corrector_Desired : constant Boolean := True;
   --  True is good here.

   Extrapolation_Desired : constant Boolean := True;

end Predictor_2;
