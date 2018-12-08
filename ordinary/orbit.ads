
--
-- package Orbit
--
-- Package defines a data structure for the 2-body problem.
-- The problem demonstrates the importance of error control 
-- with variable step size. Integrating highly elliptical 
-- orbits (near collisions between inteacting bodies) with 
-- variable step size method is orders of magnitude faster 
-- than the fixed step size method.
--
-- Pretend here that the two body system is a classical hydrogen
-- atom in two dimensions.
--
-- Unit of distance is Bohr radius.
-- Unit of energy is the ground state energy of hydrogen
-- Unit of time planck's constant (h-bar) divided by the unit of energy.
--
-- In these units (the equivalent of) Kepler's law is
--      (w**2)*(r**3) = 4
-- where w is orbital frequency, and r is the separation of the
-- two bodies. For example, if we set r = 1 (ground state of hydr.) then
-- w = 2, and the speed is w*r = 2. Consequently a circular orbit of
-- r = 1 has velocity of magnitude 2 perpindicular to its positional
-- vector.  In the demo program, velocity << 2, so that the orbit is
-- highly elliptical.

generic

   type Real is digits <>;

package Orbit is

   type Dyn_Index is range 0..3;

   type Dynamical_Variable is array(Dyn_Index) of Real;

   DynZero : constant Dynamical_Variable := (others => (+0.0));

   function F 
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable;

   -- Defines the equation to be integrated,
   -- dY/dt = F (t, Y).  Even if the equation is t or Y
   -- independent, it must be entered in this form.

   function Energy 
     (Y : in Dynamical_Variable) 
      return Real;

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

end Orbit;
