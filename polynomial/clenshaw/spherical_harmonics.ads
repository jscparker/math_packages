--
-- package Spherical_Harmonics
--
-- Actually makes just the part that varies in Theta..leaves out
-- the Exp(i*m*Phi) / Sqrt(2*Pi), which is all that's required to complete
-- the spherical harmonic.  Uses package Clenshaw to make the polynomials.
--
-- The results may be thought of as normalized versions of Associated
-- Legendre Functions, Q(k, m, X) (which in standard form are not 
-- normalized, and are not actually polynomials if m is odd).
--
-- Designed to be most efficient when k and X change a lot, and m
-- changes rarely.  A table of Normalization factors for fixed
-- m and varying k is made.  if k is fixed and m varies a lot, then
-- this method is grossly inefficient.
--
-- Package should be thought of as an example of how to use package
-- Clenshaw and package A_Legendre to make spherical harmonics.  There
-- will always be more efficient methods. Best method is problem dependent.
--
-- (un-normalized) Spherical_Harmonic (m, k):   ( l = k + m )
--
--     Q_0 (m, X)      = (-1)**m * Sqrt(1-X*X)**m
--     Q_1 (m, X)      = X * (2*(m+1) - 1) * Q_0 (m, X) = Alpha*Q_0
--     Q_k (m, X)      = X * ((2*(m+k) - 1) / k) * Q_k-1 (m, X)
--                          -((k + 2*m - 1) / k) * Q_k-2 (m, X)
--     Alpha (k, m, X) = X * (2*(m+k) - 1) / k
--     Beta (k, m, X)  = -(k + 2*m - 1) / k
--
-- Functions are orthogonal on the interval [-1,1] with
-- weight function W(X) = 1.  Orthogonality is respect integration, not
-- summation of discrete data points.  Normalizing integral:
--
--     Int (Q_k(m, X) * Q_k(m, X) * W(X))
--                 = 2*(k+2*m)! / ((2*k + 2*m + 1) * (k!*(2m-1)!!**2))
--
-- The actual Associated Legendre Polys are usually defined with (2m-1)!! times
-- the Q_0 given above, but this leads to overflow, so its put in the Norm. 
--
-- The m values for the Associated Legendre Polys are always non-negative. When
-- you use Associated Legendre Polys to make spherical, (where m is in -l..l)
-- then use Abs(m) to make the Assoc. Legendre Polys.
-- Azimuthal quantum number l is equal to k + m.
-- Notice that l_min = m, (since m_max = l).  k is in 0..inf,
-- so l = k + m is in m..inf.  So in the standard l and m notation, (m > -1):
--
--     Int (Q_k(m, X) * Q_k(m, X) * W(X))
--                 = 2*(l+m)! / ((2*l + 1) * (l-m)! * (2m-1)!!))
--
-- To make spherical harmonics, which
-- are normalized respect integration in Phi also, you multiply the results
-- of Clenshaw (the Associated Legendre Functions) by 
-- Sqrt [(2*l + 1) * (l-m)! * (2m-1)!! / 2*(l+m)!] and by
-- [1.0 / Sqrt (2*Pi)] * Exp(i*m*Psi).
--
generic

   type Real is digits <>;

   with function Sqrt (X : Real) return Real;
   with function  Exp (X : Real) return Real;
   with function  Log (X : Real) return Real;
   
   type Base_Poly_Index is range <>;

   Poly_Limit : Base_Poly_Index;
   
package Spherical_Harmonics is

   pragma Assert (Poly_Limit > 1 and Base_Poly_Index'First <= 0);
   -- Index l is of type Base_Poly_Index, and l must be >= 0.
 
   function X_Lower_Bound return Real;  -- -1.0
   function X_Upper_Bound return Real;  -- +1.0

    --  The Spherical_Harms use a table to speed up calculation: designed 
    --  to be efficient when X and k change a lot, and m changes rarely.

   function Spherical_Harm 
     (l          : Base_Poly_Index;
      m          : Real;
      X          : Real;
      Iterations : Positive := 1)      
      return Real;
   --  Uses Evaluate_Qs() to get array of Spherical_Harm(X) in k=[0, l].
   --
   --  For efficiency, ought to add another version of Spherical_Harm that
   --  returns the whole array of Vals, from l = 0 to l_max.
   --
   --  Seems to work ok to l,m >> 3000.
   --
   --  Seems to work to machine precision to very high l,m, so haven't
   --  found a limit where Iterations > 1 is needed.

   function Spherical_Harm_2
     (l          : Base_Poly_Index;
      m          : Real;
      X          : Real;
      Iterations : Positive := 1)      
      return Real;
   --  Uses Clenshaw's summation to get Spherical_Harmonics.  Used for
   --  testing, mainly. 

   function Poly_Weight (X : Real) return Real;

end Spherical_Harmonics;
