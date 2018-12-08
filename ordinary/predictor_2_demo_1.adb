
with Predictor_2;
with Text_IO; use Text_IO;
with Sinu_2;
with Ada.Numerics.Generic_Elementary_Functions;

procedure predictor_2_demo_1 is

   type Real is digits 15;

   package Sinusoid is new Sinu_2 (Real);
   use Sinusoid;

   package Sin_Integrate is new 
      Predictor_2 (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-");
   use Sin_Integrate;

   package mth is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use mth;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

   Previous_Y, Newest_Y         : Dynamical_Variable;
   Previous_Y_dot, Newest_Y_dot : Dynamical_Variable;
   Final_Time, Starting_Time : Real;
   Previous_t, Newest_t      : Real;
   Delta_t                   : Real := 1.0 / 8.0;

   Steps : Real := 10_000.0;

begin
   new_line;
   put ("Test integrates an equation whose solution is Sin(t).");
   new_line;
   put ("The equation is (d/dt)**2 (Y) = -Y. (Y = Sin(t) has a period of 2*Pi.)");
   new_line(2);
   put ("Enter number of time steps (try 512): ");
   get (Steps);
   new_line;
   put ("Every time the integration advances this number of steps, ERROR is printed.");
   new_line;

   --  choose initial conditions  

   Starting_Time        :=  0.0;
   Final_Time           :=  64.0;

   Previous_Y(0)     := 0.0;
   Previous_Y_dot(0) := 1.0;
   Previous_Y(1)     := 0.0; --unused... just too lazy to remove it.
   Previous_Y_dot(1) := 0.0; --unused
   Previous_t := Starting_Time;
   Delta_t    := Final_Time - Starting_Time;

   for i in 1..28 loop

      Newest_t := Previous_t + Delta_t;

      Integrate
        (Final_Y            => Newest_Y,      -- the result (output).
         Final_deriv_Of_Y   => Newest_y_dot,
         Final_Time         => Newest_t,      -- integrate out to here (input).
         Initial_Y          => Previous_Y,    -- input an initial condition (input).
         Initial_deriv_Of_Y => Previous_Y_dot,
         Initial_Time       => Previous_t,    -- start integrating here (input).
         No_Of_Steps        => Steps);        -- use this no_of_steps

      Previous_Y     := Newest_Y;
      Previous_Y_dot := Newest_y_dot;
      Previous_t     := Newest_t;

    --put ("Final value of numerically integrated sin:  ");
    --put (Newest_Y(0), Aft => 7);
      new_line;
      put ("Time = t =");
      put (Newest_t, Aft => 7);
      put (",     Error = Sin(t)-Y = ");
      put (Sin (Real (Newest_t)) - Newest_Y(0), Aft => 7);

   end loop;

   new_line(2);
   put ("Total elapsed time is:");
   put (Newest_t / (2.0*3.14159266), Aft => 7);
   put (" in units of 2 pi.");

   new_line(2);
   put ("Final Remark: this (20th order) Predictor-Corrector requires half as many");
   new_line;
   put ("time steps per interval as the 17th order Predictor-Corrector to get near");
   new_line;
   put ("machine precision. That means it requires about 1/12th as many evaluations");
   new_line;
   put ("of F(t, Y) as the 8th order Runge-Kutta, (for this class of equation anyway).");
   new_line;
   put ("If evaluation of F(t, Y) dominates overhead at run-time, as it does in most");
   new_line;
   put ("N-body problems, then this specialized Predictor-Corrector is a good idea");
   new_line;
   put ("here. In fact Predictor-Correctors are still commonly used in these problems.");
   new_line;
   put ("Of course the 20th order Predictor-Corrector is restricted to a special");
   new_line;
   put ("class of differential equation: (d/dt)^2 Y = F(t, Y).");
   new_line;

   new_line(2);

end;
