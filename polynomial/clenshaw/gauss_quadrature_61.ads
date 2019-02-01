
-- package Gauss_Quadrature_61
--
-- Package provides a 61 pt. Gauss-Kronrod quadrature rule.  Tables are
-- good to about 30 significant figures.
-- A 30-point Gaussian rule is used for error estimate.(From Quadpack.)
--
-- To estimate the idefinite integral of a function F(X) in an interval
--  (X_Starting, X_Final), use the following fragment:
--
--	Find_Gauss_Nodes (X_Starting, X_Final, X_gauss);
--
--	for I in Gauss_Index_61 loop
--	   F_val(I) := F (X_gauss(I));
--	end loop;
--      --  Evaluate the function at the 61 callocation points just once.
--
--      Get_Integral (F_val, X_Starting, X_Final, Area, Rough_Error);
--
-- The array X_gauss(I) is the set of points X that the
-- function F is evaluated at. The function is multiplied by
-- the Gauss_Weights associated with each of these points, and
-- result is summed to give the area under the curve in the
-- interval defined above.
--
-- Rough_Error is usually wrong by orders of magnitude.  But
-- it's usually larger than the actual error. Works best when
-- error is near machine epsilon.
--
generic

   type Real is digits <>;

   with function Sqrt(X : Real) return Real is <>;

package Gauss_Quadrature_61 is

   type Gauss_Index_61 is range -30 .. 30;

   type Gauss_Values is array (Gauss_Index_61) of Real;

   --  These are the 61 callocation points, the X axis points at
   --  which the function is evaluated.  They are calculated by
   --  the following routine:

   procedure Find_Gauss_Nodes
     (X_Starting : in     Real;
      X_Final    : in     Real;
      X_gauss    :    out Gauss_Values);

   type Function_Values is array(Gauss_Index_61) of Real;

   --  Fill array F_val of this type with values function F: 
   --
   --   for I in Gauss_Index_61 loop
   --      F_val(I) := F (X_gauss(I));
   --   end loop;

   procedure Get_Integral
     (F_val       : in     Function_Values;
      X_Starting  : in     Real;
      X_Final     : in     Real;
      Area        :    out Real;
      Rough_Error :    out Real);

   procedure Gauss_61_Coeff_Test;

end Gauss_Quadrature_61;
