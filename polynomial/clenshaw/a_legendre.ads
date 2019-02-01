
-- package A_Legendre
--
-- Data structure for Associated Legendre Polynomials. Used by package
-- Clenshaw to generate the Associated Legendre functions via recurrance 
-- relations.
--
-- The Norms of the functions are calculated separately from the functions.
-- That's so that calculation of the norms can be moved outside the inner
-- loop that generates the functions.  Calculating the Norms may be 
-- excessively expensive in time.
--
-- Function Norm is so slow it should be used to fill a table, rather 
-- than called excessively.
--
-- Up to roughly order 900(?) is OK for k and m.  Beyond that still more
-- tricks are needed, but the tricks vary with l and m, so I don't bother
-- with them.
--
-- If you enter k < 0, or m < 0, then Constraint_error is raised.  If you
-- are using this data structure to make Spherical Harmonics, (which allow
-- m < 0), then first set m = Abs (m);  the sign of m influences only the
-- Exp (i*m*Phi) part of the Spherical Harmonic.
--
-- (un-normalized) Associated Legendre (m, k):   ( l = k + m )
--
--     Q_0 (m, X) = (-1)**m * Sqrt(1-X*X)**m
--     Q_1 (m, X) = X * (2*(m+1) - 1) * Q_0 (m, X) = Alpha*Q_0
--     Q_k (m, X) = X * ((2*(m+k) - 1) / k) * Q_k-1 (m, X)
--                     -((k + 2*m - 1) / k) * Q_k-2 (m, X)
--     Alpha (k, m, X) = X * (2*(m+k) - 1) / k
--     Beta (k, m, X)  = -(k + 2*m - 1) / k
--
-- Functions are orthogonal on the interval [-1,1] with
-- weight function W(X) = 1.  Orthogonality is respect integration, not
-- summation of discrete data points.  Normalizing integral:
--
--     Integral (Q_k(m, X) * Q_k(m, X) * W(X))
--            = (k+2*m)! / ((k + m + 0.5) * k! * (2m-1)!!**2)
--
-- The actual Assoc. Legendre Functions are usually defined with (2m-1)!! times
-- the Q_0 given above, but this leads to overflow, so it's put in the Norm. 
-- The m values for the Assoc. Legendre Polys are always non-negative. When
-- you use Assoc. Legendre Polys to make spherical, (where m is in -l..l)
-- then use Abs(m) to make the Associated Legendre Functions.
--
-- Data structure for instantiation of Clenshaw:
--
generic
 
   type Real is digits <>;

   with function Sqrt (X : Real) return Real;
   with function  Exp (X : Real) return Real;
   with function  Log (X : Real) return Real;
   
   type Base_Poly_ID is range <>;
   --  Must include 0 in its range.  This is checked???
    
package A_Legendre is

   function X_Lower_Bound return Real;  -- -1.0
   function X_Upper_Bound return Real;  -- +1.0

   function Alpha (k : Base_Poly_ID; m : Real; X : Real) return Real;

   function Beta (k : Base_Poly_ID; m : Real; X : Real) return Real;

   function Q_0 (m : Real; X : Real) return Real;

   function Normalization_Factor (k : Base_Poly_ID; m : Real) return Real;
   -- Multiply the Q's by this to normalize.

   function Normalization_Factor_0 (k : Base_Poly_ID; m : Real) return Real;
   -- Alternative norm; for testing.

   function Poly_Weight (X : Real) return Real;

end A_Legendre;
