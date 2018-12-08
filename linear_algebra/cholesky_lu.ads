
-- PACKAGE Cholesky_LU
--
-- Cholesky's algorithm for LU decomposition of Real positive-definite
-- square matrices.  Positive definite matrices are a special class of
-- symmetric, non-singular matrices for which no pivioting is required
-- in LU decomposition.  LU_decompose factors the matrix A into A = L*U,
-- where U = transpose(L).  LU_Decompose writes over A with L and U.
-- U is upper triangular, and is placed in the upper triangular part of A.
--
-- The procedure does not test the matrix for positive definiteness.
-- (A good way to test for positive definiteness is to run LU_decompose 
-- on it to see if Constraint_Error is raised.)  All positive definite 
-- matrices are symmetric and their central diagonals have no zero
-- valued elements or negative valued elements.
--
-- A square (N X N) matrix with elements of generic type Real
-- is input as "A" and returned in LU form.  A must be symmetric.  This
-- isn't tested.  Instead only the lower triangle of A is read.  It
-- is assumed that the upper triangle of A is the transpose of the lower.
--
-- The LU form of A can be used to solve simultaneous linear
-- equations of the form  A*X = B.  The column vector B is input
-- into procedure Solve, and the solution is returned as X.
--

generic

   type Real is digits <>;

   type Index is range <>;
   --  Defines the maximum size that the matrix can be.
   --  The maximum matrix will be (N X N) where
   --  N = Index'Last - Index'First + 1.  This is the storage
   --  set aside for the matrix, but the user may use the routines
   --  below on matrices of any size up to and including (N X N).

   type Matrix is array(Index, Index) of Real;
   --  Row major form is appropriate for Matrix * Col_Vector_Vector
   --  operations, which dominate the algorithm in procedure Solve.

package Cholesky_LU is

  type Row_Vector is array(Index) of Real;

  subtype Col_Vector is Row_Vector;

  procedure LU_Decompose 
    (A               : in out Matrix;                  -- A is over-written with LU.
     Diag_Inverse    :    out Col_Vector;
     Final_Index     : in     Index   := Index'Last;
     Starting_Index  : in     Index   := Index'First);
      
  --  Destroys A and writes over it with the LU decomposition of A.
  --  Only operates in range Starting_Index .. Final_Index.
  --  A must be symmetric (not checked). 
  --  Constraint_Error is raised if A is not positive definite.
  --  Constraint_Error is raised if evidence of singularity is detected.
  --  In both cases this evidence may be due to numerical error.

  procedure Solve 
    (X              :    out Col_Vector;
     B              : in     Col_Vector;
     A_LU           : in     Matrix;
     Diag_Inverse   : in     Col_Vector;
     Final_Index    : in     Index := Index'Last;
     Starting_Index : in     Index := Index'First);
     
  --  Solves for X in the equation A*X = b.  You enter the LU
  --  decomposition of A, not A itself.  A_Lu and Diag_Inverse are
  --  the objects returned by LU_decompose.
 
  function Product
    (A              : in Matrix;
     X              : in Col_Vector;
     Final_Index    : in Index := Index'Last;
     Starting_Index : in Index := Index'First) 
     return Col_Vector;
     
  --  Matrix vector multiplication.
  
  function "-"(A, B : in Col_Vector) return Col_Vector;

private

  Zero : constant Real := (+0.0);
  One  : constant Real := (+1.0);
  Two  : constant Real := (+2.0);

  Min_Allowed_Real : constant Real := Two ** (Real'Machine_Emin + 4);
  --  Must be positive.  If pivots are found to be smaller than this, then
  --  it is taken as evidence that the matrix is not postive definite,
  --  and Constraint_Error is raised.

end Cholesky_LU;
