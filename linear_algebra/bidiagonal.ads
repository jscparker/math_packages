
--
-- procedure Bi_Diagonalize
--
-- Transforms matrix A into bi-diagonal matrix B,
--
--       U' * A * V = B
--
-- The diagonal and superdiagonals of B are non-zero. All other elements
-- of B are zero.
--
-- U and V are orthogonal matrices - constructed from the products of 2 x 2 
-- Givens rotation matrices. Here, U matrix has the same shape as A, and V
-- is a square matrix with the same number of columns a A. (V has the same
-- shape as A'*A = A_transpose * A.)
--
-- Can have more Rows than Columns, but not more Columns than Rows.
--

generic

   type Real is digits <>;

   type Col_Index is range <>;
   type Row_Index is range <>;
   -- Must have Col_Index'First = Row_Index'First.  This is checked.
   -- (Want the diagonal to be A(i, i).)

   type A_Matrix is array(Row_Index, Col_Index) of Real;

   type U_Matrix is array(Row_Index, Row_Index) of Real;

   type V_Matrix is array(Col_Index, Col_Index) of Real;

package Bidiagonal is

   function U_Identity return U_Matrix;

   function V_Identity return V_Matrix;

   -- Row_Index must contain Col_Index in its range:

   pragma Assert (Row_Index (Col_Index'First) >= Row_Index'First);
   pragma Assert (Row_Index (Col_Index'Last)  <= Row_Index'Last);

   -- A is the matrix to be transformed with a series of Givens rotations:
   -- A = U * B * V'.  U and V are orthogonal matrices.
   --
   -- The procedure only operates on square matrices, or matrices with more 
   -- rows than columns.
   --
   -- Procedure checks: pragma Assert (No_of_Rows >= No_of_Cols);
   --
   -- Initial_V_Matrix = Identity, unless you preprocessed A with an L*Q
   -- transform, in which case  Initial_V_Matrix = Q, and use L for argument A.
   --
   -- Starting_Row is set to Starting_Col:

   procedure Bi_Diagonalize
     (A                : in out A_Matrix;  -- A becomes the B in A = U * B * V'
      V                :    out V_Matrix;  -- Need to input Initial_V_Matrix below.
      U                :    out U_Matrix;  -- Initialized with Identity.
      Initial_V_Matrix : in     V_Matrix;  -- Normally must use function V_Identity.
      Initial_U_Matrix : in     U_Matrix;  -- Normally just use function U_Identity.
      Starting_Col     : in     Col_Index  := Col_Index'First;
      Final_Col        : in     Col_Index  := Col_Index'Last;
      Final_Row        : in     Row_Index  := Row_Index'Last;
      Matrix_U_Desired : in     Boolean    := True;
      Matrix_V_Desired : in     Boolean    := True);

   procedure Zero_Shift_Bidiagonal_QR
     (A                : in out A_Matrix;  -- A becomes the B in A = U * B * V'
      V                : in out V_Matrix;  -- 
      U                : in out U_Matrix;  -- Initialized with Identity
      Starting_Col     : in     Col_Index  := Col_Index'First;
      Final_Col        : in     Col_Index  := Col_Index'Last;
      Final_Row        : in     Row_Index  := Row_Index'Last;
      Matrix_U_Desired : in     Boolean    := True;
      Matrix_V_Desired : in     Boolean    := True);

end Bidiagonal;
