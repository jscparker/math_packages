
-- PACKAGE Runge_5th
--
-- Package implements the 5th order Runge Kutta of Prince and Dormand.
--
-- The program integrates     dY/dt = F (t, Y)    where t is the
-- independent variable (e.g. time), vector Y is the dependent 
-- variable, and F is some function. 
--
-- When high accuracy is required, higher order routines 
-- are a good idea; but if very low accuracy is fine, then the 5th order
-- (6 stage) Runge-Kutta might be better. If the equation is
-- is stiff, or if you are forced to take very small step sizes for
-- any other reason, then the lower order method is probably faster.
--
-- NOTES ON USE
--
-- (see runge_8th.ads)

generic

   type Real is digits <>;
   --  The independent variable has type Real.  Its called Time
   --  throughout the package as though the differential
   --  equation dY/dt = F (t, Y) were time dependent.

   type Dyn_Index is range <>;

   type Dynamical_Variable is array(Dyn_Index) of Real;
   --  The dependent variable.
   --
   --  To use this routine, reduce higher order differential
   --  equations to 1st order.
   --  For example a 3rd order equation in X becomes a first order
   --  equation in the vector Y = (X, dX/dt, d/dt(dX/dt)).

   with function F
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines the equation to be integrated: dY/dt = F (t, Y).
   --  Even if the equation is t or Y
   --  independent, it must be entered in this form.

   with function "*"
     (Left  : Real;
      Right : Dynamical_Variable)
      return Dynamical_Variable is <>;
   --  Defines multiplication of the independent by the
   --  dependent variable. An operation of this sort must exist by
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

package Runge_5th is

   type Step_Integer is range 1 .. 2**31-1;
   --  The step size Delta_t is never input. Instead
   --  you input Starting_Time, Final_Time, and No_Of_Steps.
   --  If there is no error control then Delta_t is
   --  (Final_Time - Starting_Time) / No_Of_Steps.  If there is
   --  error control, then the Delta_t given above is the
   --  initial choice of Delta_t, but after the first step Delta_t
   --  becomes variable.

   procedure Integrate
     (Final_State           :    out Dynamical_Variable;
      Final_Time            : in     Real;
      Initial_State         : in     Dynamical_Variable;
      Initial_Time          : in     Real;
      No_Of_Steps           : in     Step_Integer;
      Error_Control_Desired : in     Boolean       := False;
      Error_Tolerance       : in     Real          := 1.0E-6);


private

   Zero : constant := +0.0;
   Half : constant := +0.5;
   One  : constant := +1.0;
   Two  : constant := +2.0;
   Four : constant := +4.0;
   One_Sixteenth  : constant := +0.0625;
   One_Eighth  : constant := +0.125;
   One_Fifth   : constant := +0.2;
   Four_Fifths : constant := +0.8;

   Test_of_Runge_Coeffs_Desired : constant Boolean := False;
   --  Test exported by Runge_Coeffs_pd_5 prints error msg if failure detected.
   --  Make true during testing.

end Runge_5th;
