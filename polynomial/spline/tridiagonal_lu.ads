
-- package Tridiagonal_LU
--
-- The package implements Crout's method for LU decomposition of
-- tri-diagonal matrices.  Matrix A is input in the form of three
-- diagonals: the central diagonal of A, indexed by 0, and the
-- two side diagonals indexed by -1 and 1.
--
-- The LU form of A can then be used to solve simultaneous linear
-- equations of the form  A * X = B.  The column vector B is input
-- into procedure Solve, and the solution is returned as X.

generic

   type Real is digits <>;

   type Index is range <>;

package Tridiagonal_LU is

   type Diagonal is array(Index) of Real;
   D : constant := 1;
   type DiagonalID is range -D..D;
      -- The lower diagonal is the -1; The upper diagonal is the +1.

   type Matrix is array(DiagonalID) of Diagonal;
   -- Row major form is appropriate for Matrix*ColumnVector
   -- operations, which dominate the algorithm in procedure
   -- Solve.
   -- The lower diagonal is the -1; The upper diagonal is the +1.

   procedure LU_Decompose
     (A            : in out Matrix;
      Index_Start  : in     Index := Index'First;
      Index_Finish : in     Index := Index'Last);

   -- In the output matrix A, the lower triangular matrix L is stored
   -- in the lower triangular region of A, and the upper, U, is stored
   -- in the upper triangular region of A.
   -- The diagonal of U is assumed to be entirely 1.0, hence the output
   -- matrix A stores the diagonal elements of L along its diagonal.
   -- The matrix to be decomposed is (M X M) where
   -- M = Index_Finish - Index'First + 1.  The Matrix A will be much larger
   -- if Index'Last > Index_Finish, but all values of A with row or column
   -- greater than Index_Finish are ignored.

   subtype Column is Diagonal;

   procedure Solve
     (X            :    out Column;
      A            : in     Matrix;
      B            : in     Column;
      Index_Start  : in     Index := Index'First;
      Index_Finish : in     Index := Index'Last);

   -- Solve for X in the equation A X = B. The matrix A is input
   -- in LU form. Its top triangular part is U, and its lower triangular
   -- is L, where L*U = A. The diagonal elements of A hold the diagonal
   -- elements of U, not L.  The diagonal elements of L are assumed to
   -- equal 1.0. The output of LU_Decompose is in suitable form for "Solve".

   matrix_is_singular : exception;

   Epsilon : Real := 16.0 * Real'Safe_Small;

   Set_Zero_Valued_Pivots_To_Epsilon : constant Boolean := False;
   -- If set to true then the pivot is given the value
   -- epsilon, and no exception is raised when matrix is singular.

end Tridiagonal_LU;
