--
-- PACKAGE Chebychev_Quadrature
--
-- The package provides tables and procedures for Gaussian Quadrature
-- based on Chebychev polynomials (Gaussian-Chebychev Quadrature).
-- The number of points used in the quadrature (callocation points)
-- is the generic formal parameter: No_Of_Gauss_Points.
--
-- Usually not the best quadrature routine, but works well in special
-- cases.
--
-- Use the following fragment to calculate the definite integral
-- of a function F(X) in an interval (X_Starting, X_Final):
--
--	Find_Gauss_Nodes (X_Starting, X_Final, X_gauss);
--
--	for I in Gauss_Index loop
--	   F_val(I) := F (X_gauss(I));
--	end loop;
--      --  Evaluate the function at the callocation points just once.
--
--      Get_Integral (F_val, X_Starting, X_Final, Area);
--
-- Notice that this package can be instantiated with extended precision
-- floating point - recommended if machine precision is necessary.
-- (The extended precision weights and Q-points are rounded to
-- normal precison after calculation, but its not terribly important.)
--
generic

   type Real is digits <>;

   No_Of_Gauss_Points : Positive;
   
   with function Cos (X : Real) return Real is <>;
   with function Sin (X : Real) return Real is <>;

package Chebychev_Quadrature is

   type Gauss_Index is new Positive range 1 .. No_Of_Gauss_Points;
 
   type Gauss_Values is array(Gauss_Index) of Real;

   -- Plug X_Starting and X_Final into the following to get the
   -- the points X_gauss at which the function to
   -- be integrated must be evaluated:
 
   procedure Find_Gauss_Nodes
     (X_Starting : in     Real;
      X_Final    : in     Real;
      X_gauss    :    out Gauss_Values);

   type Function_Values is array(Gauss_Index) of Real;

   -- Fill array F_val of this type with values function F: 
   --
   --   for I in Gauss_Index loop
   --      F_val(I) := F (X_gauss(I));
   --   end loop;

   procedure Get_Integral
     (F_val       : in     Function_Values;
      X_Starting  : in     Real;
      X_Final     : in     Real;
      Area        :    out Real);

end Chebychev_Quadrature;

