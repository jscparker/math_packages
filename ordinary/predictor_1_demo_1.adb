
with Predictor_1;
with Text_IO; use Text_IO;
with Sinu;
with Ada.Numerics.Generic_Elementary_Functions;

procedure predictor_1_demo_1 is

   type Real is digits 15;

   package Sinusoid is new Sinu (Real);
   use Sinusoid;

   package Sin_Integrate is new 
      Predictor_1 (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-", Norm);
   use Sin_Integrate;

   package math is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use math;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Step_Integer);
   use iio;

   Initial_Condition         : Dynamical_Variable;
   Previous_Y, Newest_Y      : Dynamical_Variable;
   Final_Time, Starting_Time : Real;
   Previous_t, Newest_t      : Real;
   Delta_t                   : Real;

   Steps : Step_Integer;

begin
   new_line;
   put ("Test integrates an equation whose solution is Sin(t).");
   new_line;
   put ("The equation is (d/dt)**2 (Y) = -Y. (Y = Sin(t) has a period of 2*Pi.)");
   new_line(2);
   put ("Enter number of time steps (try 1024): ");
   get (Steps);
   new_line;
   put ("Every time the integration advances this number of steps, ERROR is printed.");
   new_line;

   --  choose initial conditions  

   Initial_Condition(0) :=  0.0;
   Initial_Condition(1) :=  1.0;

   Starting_Time        :=  0.0;
   Final_Time           :=  64.0;

   Previous_Y := Initial_Condition;
   Previous_t := Starting_Time;
   Delta_t    := Final_Time - Starting_Time;

   for i in 1..28 loop

      Newest_t := Previous_t + Delta_t;

      Integrate
         (Final_State        => Newest_Y,   -- the result (output).
          Final_Time         => Newest_t,   -- integrate out to here (input).
          Initial_State      => Previous_Y, -- input an initial condition (input).
          Initial_Time       => Previous_t, -- start integrating here (input).
          No_Of_Steps        => Steps);     -- use this no_of_steps

      Previous_Y := Newest_Y;
      Previous_t := Newest_t;

    --put ("Final value of numerically integrated sin:  ");
    --put (Newest_Y(0), Aft => 7);
      new_line;
      put ("Time = t =");
      put (Newest_t, Aft => 7);
      put (",     Error = Sin(t)-Y = ");
      put (Sin (Newest_t) - Newest_Y(0), Aft => 7);

   end loop;

   new_line(2);
   put ("Total elapsed time is:");
   put (Newest_t / (2.0*3.14159266), Aft => 7);
   put (" in units of 2 pi.");

   new_line(2);
   put ("Final Remark: the Predictor-Corrector does 2 evaluations of F(t, Y)");
   new_line;
   put ("each time step. The Runge-Kutta does 13. So the Predictor-Corrector"); 
   new_line;
   put ("is the more efficient if F(t, Y) is slow to evaluate. The gravitational");
   new_line;
   put ("N-body problem (large N) is a classic example.  On the other hand, if you");
   new_line;
   put ("want to use a large time step (and if the larger error is unimportant)");
   new_line;
   put ("then the Runge-Kutta may be better. The stability of the Predictor-");
   new_line;
   put ("Corrector is so poor that it restricts you to small time steps.");

   new_line(2);

end;
