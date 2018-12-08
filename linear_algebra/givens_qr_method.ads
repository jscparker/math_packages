generic

   type Real is digits <>;

   type Index is range <>; 

   type A_Matrix is array(Index, Index) of Real;

package Givens_QR_Method is

   subtype C_Index is Index;
   subtype R_Index is Index;

   function Identity return A_Matrix;

   -- Input matrix A (in Upper Hessenberg form) and also matrix Q.
   --
   -- The original A was transformed by: 
   --
   --    A_original = Q * A * Q_transpose    where  A  is hessenberg.
   --
   -- (i.e., want eventually: Q_transpose * A_original * Q = A = A_diagonal, so
   -- that the column vectors of Q (called Z in Peters_Eigen) are eigvecs of A.)
   --
   -- Let G be a Givens rotation matrix. (Actually it will be a series of them.)
   --
   -- We develop the decomposition further by inserting G_transpose * G and
   -- its transpose on both sides of A in Q * A * Q_transpose:
   --
   --   (Q_transpose * G_transpose) * (G * A * G_transpose) * (G * Q)
   --
   --     =  Q_new_transpose * A_new * Q_new
   --
   -- So to develop A we find
   --
   --   A := G * A * G_transpose
   --
   -- And to develop Q we find
   --
   --   Q := G * Q
   --
   -- (With shift "s", shift A first: B = A - s * I;  
   --  Use QR to calculate G*R where R is upper triangular:
   --  
   --   G * R = B
   --   R * G = G' * B * G
   --   R * G + s * I = G' * A * G
   --  
   -- So the matrix returned as the new A is R*G + s*I.)

   -- works only for Upper Hessenberg matrices A:

   procedure Lower_Diagonal_QR_Iteration
     (A            : in out A_Matrix;
      Q            : in out A_Matrix;
      Shift        : in     Real;
      Starting_Col : in     C_Index  := C_Index'First;
      Final_Col    : in     C_Index  := C_Index'Last);

end Givens_QR_Method;
