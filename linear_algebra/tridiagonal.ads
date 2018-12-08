
---------------------------------------------------------------------------
-- package Tridiagonal, symmetric matrix tridiagonalization
-- Copyright (C) 2018 Jonathan S. Parker
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
---------------------------------------------------------------------------

-- package Tridiagonal
--
-- Tridiagonalization of symmetric matrices, using Givens rotations.

generic

   type Real is digits <>;

   type Index is range <>; 

   type A_Matrix is array(Index, Index) of Real;

package Tridiagonal is

   subtype C_Index is Index;
   subtype R_Index is Index;

   function Identity return A_Matrix;

   -- The symmetric input matrix A is transformed with similarity transformations:
   --
   --      A_tridiagonal = Q_transpose * A * Q.
   --
   -- The Q's are orthogonal matrices constructed from the products of
   -- 2 x 2 Givens matrices.
   -- Q matrix has the same shape as A.

   procedure Tridiagonalize
     (A                : in out A_Matrix;
      Q                :    out A_Matrix;
      Starting_Col     : in     C_Index  := C_Index'First;
      Final_Col        : in     C_Index  := C_Index'Last;
      Initial_Q        : in     A_Matrix := Identity;
      Q_Matrix_Desired : in     Boolean  := True);

   --  Works only on Tridiagonal matrices, and returns silently if Mat_Size<3:

   procedure Lower_Diagonal_QR_Iteration
     (A                : in out A_Matrix;
      Q                : in out A_Matrix;
      Shift            : in     Real;
      Final_Shift_Col  : in     C_Index  := C_Index'Last;
      Starting_Col     : in     C_Index  := C_Index'First;
      Final_Col        : in     C_Index  := C_Index'Last;
      Q_Matrix_Desired : in     Boolean  := True);

end Tridiagonal;
