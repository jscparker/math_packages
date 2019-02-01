
-- PACKAGE Crout_LU
--
-- LU decomposition and linear equation solving, with partial pivoting,
-- for square, real valued matrices.
--
-- A square (N X N) matrix with elements of generic type Real
-- is input as "A" and returned in LU form, along with an array
-- containing information on the permutation of the rows of the
-- matrix that occurred during pivoting.
--
-- The decomposition can be performed on arbitrary diagonal blocks of A.
--
-- If Scaling_Desired = True, then matrices are scaled prior to LU
-- decomposition. First all the columns are scaled to near unity (in
-- 1-norm). Then the process is repeated for the rows.  It's not fast
-- but occasionally improves the decomposition.
--
-- The LU form of A can then be used to solve simultaneous linear
-- equations of the form  A X = B.  Column vector B is input
-- into procedure Solve, and the solution is returned as X.
--

generic

   type Real is digits <>;

   type Index is range <>;
   --  Defines the maximum size of the matrix: (N X N) where
   --  N = Index'Last - Index'First + 1.  This is storage
   --  set aside for the matrix, but the routines operate on
   --  arbitrary diagonal blocks of the Matrix.

   type Matrix is array(Index, Index) of Real;

package Crout_LU is

   type Row_Vector is array(Index) of Real;

   subtype Col_Vector is Row_Vector;

   type Rearrangement is array(Index) of Index;

   type Scale_id is (Diag_Inverse, For_Rows, For_Cols);

   type Scale_Vectors is array (Scale_id) of Row_Vector;

   procedure LU_Decompose
     (A                : in out Matrix;    -- A is overwritten with L and U
      Scalings         :    out Scale_Vectors;
      Row_Permutation  :    out Rearrangement;
      Final_Index      : in     Index   := Index'Last;
      Starting_Index   : in     Index   := Index'First;
      Scaling_Desired  : in     Boolean := False);
 
   -- In the output matrix A, the lower triangular matrix L is stored
   -- in the lower triangular region of A, and the upper, U, is stored
   -- in the upper triangular region of A.
   -- The diagonal of L is assumed to be entirely 1.0, so the output
   -- matrix A stores the diagonal elements of U along its diagonal.
   -- The matrix to be decomposed is (M X M) where
   -- M = Final_Index - Index'First + 1.  The Matrix A can be much larger
   -- than M X M, but all values of A with row or column
   -- greater than Final_Index are ignored.


   procedure LU_Solve
     (X                :    out Row_Vector;
      B                : in     Row_Vector;
      A_LU             : in     Matrix;
      Scalings         : in     Scale_Vectors;
      Row_Permutation  : in     Rearrangement;
      Final_Index      : in     Index := Index'Last;
      Starting_Index   : in     Index := Index'First);
 
   -- Solve for X in the equation A X = B. The matrix LU_of_A is the LU decomp
   -- of matrix A. Its top triangular part is U, and its lower triangular
   -- is L, where L*U = A.
   -- The output of LU_Decompose is in suitable form for "Solve".


   function Product
     (A              : Matrix;
      V              : Row_Vector;
      Final_Index    : Index := Index'Last;
      Starting_Index : Index := Index'First)
      return Row_Vector;

   function "-"(A, B : in Col_Vector) return Col_Vector;

   procedure Scale_Cols_Then_Rows
     (A              : in out Matrix;
      Scalings       :    out Scale_Vectors;
      Final_Index    : in     Index := Index'Last;
      Starting_Index : in     Index := Index'First);
   --  Returns the matrix A in its scaled form.  The only purpose
   --  is for testing.

end Crout_LU;
