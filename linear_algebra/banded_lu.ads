
---------------------------------------------------------------------------------
-- package Banded_LU, LU decomposition, equation solving for banded matrices
-- Copyright (C) 1995-2018 Jonathan S. Parker 
-- 
-- Permission to use, copy, modify, and/or distribute this software for any
-- purpose with or without fee is hereby granted, provided that the above
-- copyright notice and this permission notice appear in all copies.
-- THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
-- WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
-- MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
-- ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
-- WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
-- ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
-- OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
---------------------------------------------------------------------------------

-- PACKAGE Banded_LU
--
-- LU decomposition and linear equation solving for banded matrices. 
-- Uses Crout's method for LU decomposition of real valued matrices.
-- There's no pivoting, so iterative refinement is provided.
-- Iterative refinement can improve accuracy of equation solving
-- in cases in which the matrices are not diagonally dominant.
--
-- The decomposition can be performed on diagonal sub-blocks of
-- the matrix. (The block size must be > bandwidwidth of the matrix
-- or Constraint_Error is raised.)
--
-- Diagonally dominant matrices are a special class that doesn't need
-- pivoting in LU decomposition.  A matrix is diagonally dominant if and
-- only if the absolute val of each diagonal element is > or = the sum
-- of the absolute vals of each off-diagonal on the same row as that diagonal
-- element.
--
-- In many cases non diagonally-dominant matrices are successfully
-- handled by this package, but more error accumulates than would
-- be present if we could do pivoting.
-- Procedure Refine_Solution does Newton iterations to reduce the error
-- in these solutions.  Examples are provided in the test procedure.
--
-- The banded matrix is constructed as an array of row vectors.
-- Typically, the length of the row is very much smaller than the number
-- of rows.  For example, if the matrix is tri-diagonal, then each row vector
-- has only three element.  The middle element is the central diagonal of the
-- matrix.
--
-- The LU form of A can then be used to solve simultaneous linear
-- equations of the form  A X = B.  The column vector B is input
-- into procedure Solve, and the solution is returned as X.
--
--
-- PROCEDURE LU_decompose
--
-- matrix A:
-- To correctly construct A remember that the 1st index identifies the
-- row of the matrix element.  The 2nd index of A is
-- is the Diagonal_ID, which is 0 for the central diagonal, and
-- negative for the lower diagonals. Diagonal_ID - Col - Row. 
--
-- Final_Index:
-- The routine operates on a subset of full matrix.  The subset
-- matrix starts at Starting_Index, and ends at Final_Index.
-- The matrix to be decomposed is logically (M X M) where
-- M = Final_Index - Starting_Index + 1.  It is assumed that all
-- off-diagonals beyond No_Of_Off_Diagonals are 0.0.
--
--
-- PROCEDURE Solve
--
-- Solves for X in the equation A X = B. The matrix A_LU is
-- LU of A. In other words  L*U = A, L is stored in the lower triangular
-- regions of A_LU, and U in the upper.  (The original matrix A was written
-- over with LU.)  The output of LU_Decompose is in suitable form for "Solve".
-- 
--
-- BANDED MATRIX FORMAT
--
--    You translate (Row, Col) to (I, Diagonal_id)
--    using the formula:   I = Row, and  Diagonal_id = Col - Row.
--
-- The banded matrix looks like:
--
--      Matrix (Row)(Diagonal_id)
--
-- where Diagonal_id is in 
--
--  type Diag_Index is range -No_Of_Off_Diagonals..No_Of_Off_Diagonals;
--
-- and the matrix element at (Row, Col) has Diagonal_id = Col - Row.
--
-- Here is matrix-vector multiplication:
--
--   for Row in Index loop
--      Sum    := 0.0;
--      Start  := Max(Diag_Index'First + Row, Index'First)
--      Finish := Min(Diag_Index'Last  + Row, Index'Last)
--      for Col in Start..Finish loop
--         Sum := Sum + Matrix(Row)(Col - Row) * Vector(Col);
--      end loop;
--      Result(Row) = Sum;
--   end loop;
--
-- Notice the outer index varies fastest this way, as Ada prefers.
--

generic

   type Real is digits <>;
    
   Max_Size_Of_Matrix  : Positive;
   No_Of_Off_Diagonals : Positive;
   --  If this is 1, then if the matrix is tri-diagonal.  Must be less than
   --  Max_Size_Of_Matrix.  If it equals or excedes Max_Size_Of_Matrix, then 
   --  an assertion will detect failure.  Notice 0 off_diagonals not allowed.

package Banded_LU is

   pragma Assert (No_Of_Off_Diagonals < Max_Size_Of_Matrix);
   -- Number of off-diagonals must be less than size of matrix.

   subtype Index      is Integer range 0..Max_Size_Of_Matrix-1;
   subtype Diag_Index is Integer range 
                                 -No_Of_Off_Diagonals..No_Of_Off_Diagonals;

   --  There's arithmetic between these 2 indices, so we make them both type
   --  Integer.  We also use fact that Diagonal_ID is symmetric about 0. 

   type Row_Vector    is array(Diag_Index) of Real;
   type Banded_Matrix is array(Index)      of Row_Vector;
   --  The (Row, Col) element of the matrix is found at Matrix(Row)(Col - Row).

   type Column is array(Index) of Real;
   --  A column vector.

   procedure LU_Decompose
     (A              : in out Banded_Matrix;
      Diag_Inverse   :    out Column;
      Final_Index    : in     Index := Index'Last;
      Starting_Index : in     Index := Index'First);
 
   procedure Solve
     (X              :    out Column;
      B              : in     Column;
      A_LU           : in     Banded_Matrix;
      Diag_Inverse   : in     Column;
      Final_Index    : in     Index := Index'Last;
      Starting_Index : in     Index := Index'First);
 
   --  Solves for X in A*X = B.  A_LU is the LU decomp of A, as returned by
   --  LU_Decompose in matrix A.

   procedure Refine_Solution
     (X                :    out Column;
      B                : in     Column;
      A_LU             : in     Banded_Matrix;
      Diag_Inverse     : in     Column;
      A                : in     Banded_Matrix;
      No_Of_Iterations : in     Natural := 1;
      Final_Index      : in     Index   := Index'Last;
      Starting_Index   : in     Index   := Index'First);
    
   --  Attempts to improve the numerical quality of solution of the linear
   --  equations by use of Newton iteration. 
   --  (No pivoting is done in this version.)

   function Matrix_Val
     (A   : Banded_Matrix;
      Row : Index;
      Col : Index) 
      return Real;

   --  Given (Row, Column) the function returns the matrix element at
   --  that location.  Translates (Row, Col) to (I, Diagonal_id) using
   --  the formula I = Row, and Diagonal_id = Col - Row.
   --
   --  Banded Matrices are by definition 0 everywhere except on the diagonal
   --  bands.  So 0 is returned if (Row, Col) is not in the banded region.
 
   function Product
     (A              : in Banded_Matrix;
      X              : in Column;
      Final_Index    : in Index  := Index'Last;
      Starting_Index : in Index  := Index'First) 
      return Column;
 
   --  Multiplies a banded matrix by a column vector X.
   --  Operates only on the diagonal block bounded by 
   --      Starting_Row = Starting_Index, and Final_Row = Final_Index.

private

   Zero  : constant Real := +0.0;
   One   : constant Real := +1.0;
   Two   : constant Real := +2.0;

   Min_Allowed_Real : constant Real := Two ** (Real'Machine_Emin/4 + 2);

end Banded_LU;
