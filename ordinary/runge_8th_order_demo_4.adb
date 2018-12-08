--
-- N-Body integrations demonstrate the need for high accuracy and 
-- reliable error estimates.
--
-- N-Body equation:
--  Acceleration(Body_j)
--   = SumOver(Body_k) { Mass(Body_k) * G * DeltaR / NORM(DeltaR)**3 }
--
--   where DeltaR = (X(Body_k) - X(Body_j)) is a vector,
--   and where G is the gravitaional constant.
--
-- Natural units:
-- Unit of time is earth_orbital_period                 : t_0.
-- Unit of distance is earth's orbital semi-major axis  : a_0.
-- Unit of mass is sum of Sun and Earth mass            : m_0.
-- These are related to each other by Kepler's law:
--
--        t_0 = 2 * Pi * SQRT (a_0**3 / (G * m_0))
--
-- where G is the gravitaional constant.  If you divide variable of
-- time in the N_Body equation by t_0 to get natural units, the equation
-- becomes:
--
--  Acceleration(Body_j) = 2*2*Pi*Pi*
--          SumOver(Body_k) { Mass(Body_k) * DeltaR / NORM(DeltaR)**3 }
--
--   where DeltaR = (X(Body_k) - X(Body_j)).
--
--  Mass, time, and distance are in the natural units given above -
--  gets rid of the gravitational constant G.  In the Equation
--  as implemented below, the masses are multiplied by 2*2*pi*pi
--  when they are defined in order to avoid an excessive number of
--  multiplications during evaluation of F(t, Y).

with Runge_8th;
with Text_IO; use Text_IO;
with Four_Body;
procedure runge_8th_order_demo_4 is

   type Real is digits 15;

   package rio is new Float_IO(Real);
   use rio;
   package Bodies_4 is new Four_Body (Real);
   use Bodies_4;
   package N_Body_Solve is new 
      Runge_8th (Real, State_Index, Dynamical_Variable, F, "*", "+", "-", Norm);
   use N_Body_Solve;

   Initial_State  : Dynamical_Variable;

   Final_Y    : Dynamical_Variable;
   Previous_Y : Dynamical_Variable;
   Final_t    : Real;
   Previous_t : Real;
   ErrorTolerance : constant Real := 2.0E-12;

   -- 2nd body :

   a                   : constant Real := 21.71534093275925;
   Orbital_Period      : constant Real := 80.0;

   Orbits_Per_int      : constant Real := 10.0;
   No_Of_Steps_Per_int : constant step_integer := 120_000;
   Delta_t             : constant Real := Orbital_Period * Orbits_Per_int;

   X2, X0, Z2, Z0 : Real;

   -- Choose initial conditions. 
   -- Body1 is the larger, mass = 1.0 (sun mass).
   -- Body2 is the smaller, mass = 0.6 (sun mass).
   -- Assume the orbital period of the 2 stars is 80.0 years.
   -- Say the planet is body 3: one earth distance
   -- from the larger star with 1 earth mass.
   -- From these calculate the semimajor axis "a" in
   -- Kepler's law given above; then put the three bodies in
   -- near circular orbits and observe stability.
   -- Remember to use the reduced-mass formulas to get distance
   -- from center of mass:
   --
   -- Body_1_Radius  r1 = a*m2 / (m1+m2)
   -- Body_2_Radius  r2 = a*m1 / (m1+m2)


   --Planet_Orbital_Radius : constant Real := 1.0;-- earth's orbital period
   --Planet_Period         : constant Real := 1.0;

   --Planet_Orbital_Radius : constant Real := 2.0;-- Twice earth's orbit
   --Planet_Period         : constant Real := 2.82842712474619;

   --Planet_Orbital_Radius : constant Real := 3.0;-- Thrice earth's orbit
   --Planet_Period         : constant Real := 5.196152422706632;

   --Planet_Orbital_Radius2 : constant Real := 1.587401052; -- 1.59 earth's orbit
   --Planet_Period2         : constant Real := 2.0;-- 1.59**1.5 from Kepler's law

   Planet_Orbital_Radius1 : constant Real := 1.0;-- 1 times earth's orbit
   Planet_Period1         : constant Real := 1.0;-- 4**1.5 from Kepler's law

   Planet_Orbital_Radius2 : constant Real := 4.0; -- 1.59 earth's orbit
   Planet_Period2         : constant Real := 8.0;-- 4.0**1.5 from Kepler's law

   m1 : constant Real := Mass(0);
   m2 : constant Real := Mass(1);

   Body_1_Radius : constant Real := a * m2 / (m1+m2);
   Body_2_Radius : constant Real := a * m1 / (m1+m2);
   Body_3_Radius : constant Real := Planet_Orbital_Radius1 + Body_1_Radius;
   Body_4_Radius : constant Real := Planet_Orbital_Radius2 + Body_1_Radius;

   Body_1_Speed  : constant Real := TwoPii*Body_1_Radius / Orbital_Period;
   Body_2_Speed  : constant Real := TwoPii*Body_2_Radius / Orbital_Period;
   Ratio         : Constant Real := Body_3_Radius / Body_1_Radius;
   Body_3_Speed  : constant Real
      := (TwoPii * Planet_Orbital_Radius1 / Planet_Period1) + Body_1_Speed;
   Body_4_Speed   : constant Real
      := (TwoPii * Planet_Orbital_Radius2 / Planet_Period2) + Body_1_Speed;

   d_X, d_Z : Real;
   Orbit    : Real   := 0.0;
   Radius2  : Real;
   Min      : Real := 100000.0;
   Max      : Real := 0.0;
   Min3, Min1  : Real := 100000.0;
   Max3, Max1  : Real := 0.0;

begin

   Update_State (Initial_State, 0,
                0.0,  Body_1_Radius, Body_1_Speed,  0.0);
   Update_State (Initial_State, 1,
                0.0,   -Body_2_Radius,  -Body_2_Speed,  0.0);
   Update_State (Initial_State, 2,
                0.0,  Body_3_Radius, Body_3_Speed,  0.0);
   Update_State (Initial_State, 3,
                0.0,  Body_4_Radius, Body_4_Speed,  0.0);

   --Initial_State(0) := (0.0,  Body_1_Radius, -- 1st body XY position
   --                          Body_1_Speed,  0.0);-- 1st body UW velocity
   --Initial_State(1) := (0.0, -Body_2_Radius, -- 2nd body XY position
   --                          -Body_2_Speed, 0.0);-- 2nd body UW velocity
   --Initial_State(2) := (0.0,  Body_3_Radius, -- 3rd body XY position
   --                         Body_3_Speed,   0.0);-- 3rd body UW velocity
   --Initial_State(3) := (0.0,  Body_4_Radius, -- 4rth body XY position
   --                         Body_4_Speed,   0.0);-- 4rth body UW velocity
   -- Notice that they are all rotating clockwise
   -- looking down on the XZ plane.
 
   Previous_t := 0.0;
   Previous_Y := Initial_State;
   Final_t    := Delta_t;
 
   new_line;
 
   loop

      Integrate
        (Final_State        => Final_Y,   -- output the result.
         Final_Time         => Final_t,   -- integrate out to here.
         Initial_State      => Previous_Y, -- input an initial condition.
         Initial_Time       => Previous_t, -- start integrating here.
         No_Of_Steps        => No_Of_Steps_Per_int);

      Previous_t := Final_t;
      Previous_Y := Final_Y;
      Final_t    := Previous_t + Delta_t;


      Orbit := Orbit + Orbits_Per_int;


      X2      := State_Val (Final_Y,2,0);
      X0      := State_Val (Final_Y,0,0);
      Z2      := State_Val (Final_Y,2,1);
      Z0      := State_Val (Final_Y,0,1);
      d_X     := X2 - X0;
      d_Z     := Z2 - Z0;
      Radius2 := d_X**2 + d_Z**2;

      if Radius2 > Max then
          Max := Radius2;
      end if;
      if Radius2 < Min then
          Min := Radius2;
      end if;

      put ("Year =    "); put (Orbit * Orbital_Period); new_line;
      put ("    Max = "); put (Max); put("    Min = "); put (Min); new_line;

      X2      := State_Val (Final_Y,3,0);
      X0      := State_Val (Final_Y,0,0);
      Z2      := State_Val (Final_Y,3,1);
      Z0      := State_Val (Final_Y,0,1);
      d_X     := X2 - X0;
      d_Z     := Z2 - Z0;
      Radius2 := d_X**2 + d_Z**2;

      if Radius2 > Max3 then
          Max3 := Radius2;
      end if;
      if Radius2 < Min3 then
          Min3 := Radius2;
      end if;

      put ("    Max3 ="); put (Max3); 
      put ("    Min3 ="); put (Min3); new_line;

      X2      := State_Val (Final_Y,1,0);
      X0      := State_Val (Final_Y,0,0);
      Z2      := State_Val (Final_Y,1,1);
      Z0      := State_Val (Final_Y,0,1);
      d_X     := X2 - X0;
      d_Z     := Z2 - Z0;
      Radius2 := d_X**2 + d_Z**2;

      if Radius2 > Max1 then
          Max1 := Radius2;
      end if;
      if Radius2 < Min1 then
          Min1 := Radius2;
      end if;

      put ("    Max1 ="); put (Max1); 
      put ("    Min1 ="); put (Min1); new_line;
      --exit when Count = 2;

   end loop;

end;
