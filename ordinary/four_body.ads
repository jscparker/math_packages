--
-- package Four_Body
--
-- The package defines a data structure for the N-Body
-- gravitational problem in two dimensions.
--
-- N = 4, but can be set to anything.
--
-- N-Body equation:
--  Acceleration(Body_j)
--   = SumOver(Body_k) { Mass(Body_k) * G * Delta_R / NORM(Delta_R)**3 }
--
--   where Delta_R = (X(Body_k) - X(Body_j)) is a vector,
--   and where G is the gravitaional constant.
--
-- Natural units:
-- Unit of time is earth_orbital_period                 : t_0.
-- Unit of distance is earth's orbital semi-major axis  : a_0.
-- Unit of mass is sum of Sun and Earth mass            : m_0.
-- The 3 quantities are related to each other by Kepler's law:
--
--        t_0 = 2 * Pi * Sqrt (a_0**3 / (G * m_0))
--
-- where G is the gravitaional constant.  If you divide variable of
-- time in the N_Body equation by t_0 to get natural units, the equation
-- becomes:
--
--  Acceleration(Body_j) = 2*2*Pi*Pi*
--        SumOver(Body_k) { Mass(Body_k) * Delta_R / NORM(Delta_R)**3 }
--
--   where Delta_R = (X(Body_k) - X(Body_j)).
--
--  Mass, time, and distance are in the natural units given above,
--  which gets rid of the gravitational constant G.  In the Equation
--  as implemented below, the masses are multiplied by 2*2*pi*pi
--  when they are defined in order to avoid excessive multiplication 
--  during evaluation of F(t, Y).
--
--  Body1 is the larger, mass = 1.0.
--  Body2 is the smaller, mass = 0.6.
--  Assume that the orbital period of the 2 stars is 80.0 years.
--  Say the planet is body 3: one earth distance
--  from the larger star, with 1 earth mass.
--  From these we get the semimajor axis "a" in
--  Kepler's law given above ... then put the three bodies in
--  near circular orbits and observe stability, and
--  remember to use the reduced-mass formulas to get distance
--  from center of mass:
--
--  First_Body_Radius  r1 = a*m2 / (m1+m2)
--  Second_Body_Radius r2 = a*m1 / (m1+m2)
--
--  Planet_Orbital_Radius : constant Real := 1.0;-- earth's orbital period
--  Planet_Period         : constant Real := 1.0;
--
--  Planet_Orbital_Radius : constant Real := 2.0;-- Twice earth's orbit
--  Planet_Period         : constant Real := 2.82842712474619;
--
--  Planet_Orbital_Radius : constant Real := 3.0;-- Thrice earth's orbit
--  Planet_Period         : constant Real := 5.196152422706632;
--
-- Planet_Orbital_Radius  : constant Real := 4.0;-- 4 times earth's orbit
-- Planet_Period          : constant Real := 8.0;-- 4**1.5 from Kepler's law
--
-- a                      : constant Real := 21.71534093275925;
-- OrbitalPeriod          : constant Real := 80.0;
--

generic

   type Real is digits <>;

package Four_Body is

   No_Of_Bodies : constant := 4;

   subtype Bodies is Integer range 0 .. No_Of_Bodies-1;

   subtype XYUV_Index is Integer range 0 .. 3;

   type CartesianCoordinate is array(XYUV_Index) of Real;

   subtype State_Index is Integer range 0 .. 4*No_Of_Bodies-1;

   type Dynamical_Variable is array(State_Index) of Real;

   function F
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable;
   -- This defines the equation to be integrated as
   -- dY/dt = F (t, Y).  Even if the equation is t or Y
   -- independent, it must be entered in this form.

   procedure Update_State 
     (Y          : in out Dynamical_Variable;
      Body_id    : in     Bodies;
      X, Z, U, W : in     Real);

   function State_Val 
     (Y       : Dynamical_Variable;
      Body_id : Bodies;
      XYUV_id : XYUV_Index) return Real;

   function "*" 
     (Left  : Real;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function "+" 
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function "-" 
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function Norm (Y : Dynamical_Variable) return Real;

   pragma Inline (F, "*", "+", "-", Norm);


   R_Pi : constant := 3.14159_26535_89793_23846_26433_83279_50288;

   Pii      : constant Real := R_Pi;
   TwoPii   : constant Real := 2.0 * R_Pi;
   TwoPii2  : constant Real := 4.0 * R_Pi * R_Pi;
 
   Mass : constant array(Bodies) of Real
           := (1.0*TwoPii2, 0.6*TwoPii2, 3.0E-6 * TwoPii2, 1.0E-22 * TwoPii2);
   --  body1 has sun mass, body2 has 0.6 sun mass, and
   --  body3 has earth mass, body4 has negligible mass.

end Four_Body;
