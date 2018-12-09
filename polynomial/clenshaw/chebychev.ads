
-- package Chebychev
--
-- Data structure for generating Chebychev polynomials of the 2nd kind,
-- using Clenshaw's method:
--        Q_0(X) = 1
--        Q_1(X) = 2*X
--        Q_k(X) = 2*X*Q_k-1(X) - Q_k-2(X)
--        Alpha(k,X) =   2*X
--        Beta(k,X)  =   -1
-- They are orthogonal on the interval (-1,1) with
-- weight function W(X) = Sqrt( 1 - X*X).
-- Normalizing integral: Integral(Q_n(X) * Q_n(X) * W(X)) = Pii/2
--
generic

   type Real is digits <>;

   with function Exp (X : Real) return real;
   with function Log (X : Real) return real;
 
package Chebychev is

   X_Lower_Bound : constant Real := -1.0;
   X_Upper_Bound : constant Real := +1.0;

   -- parameter m is unused.

   type Poly_ID_Integer is new Integer;

   function Alpha (k : Poly_ID_Integer; m : Real := 0.0; X : Real) return Real;

   function Beta  (k : Poly_ID_Integer; m : Real := 0.0; X : Real) return Real;

   function Q_0 (X : Real; m : Real := 0.0) return Real;

   function Q_1 (X : Real; m : Real := 0.0) return Real;

   function Norm (k : Poly_ID_Integer; m : Real := 0.0) return Real;

   function Poly_Weight (X : Real) return Real;

end Chebychev;
