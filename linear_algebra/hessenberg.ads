
---------------------------------------------------------------------------
-- package Hessenberg
-- Copyright (C) 2011-2018 Jonathan S. Parker
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

generic

   type Real is digits <>;

   type Index is range <>; 

   type A_Matrix is array(Index, Index) of Real;

package Hessenberg is

   subtype C_Index is Index;
   subtype R_Index is Index;

   function Identity return A_Matrix;

   -- The input matrix A is transformed with similarity transformations:
   --
   --      A_hessenberg = Q_transpose * A * Q.
   --
   -- The Q's are orthogonal matrices constructed from the products of
   -- 2 x 2 Givens matrices.
   -- Q matrix has the same shape as A.

   procedure Lower_Hessenberg
     (A            : in out A_Matrix;
      Q            :    out A_Matrix;
      Starting_Col : in     C_Index := C_Index'First;
      Final_Col    : in     C_Index := C_Index'Last;
      Initial_Q    : in     A_Matrix := Identity);

   procedure Upper_Hessenberg
     (A            : in out A_Matrix;
      Q            :    out A_Matrix;
      Starting_Col : in     C_Index  := C_Index'First;
      Final_Col    : in     C_Index  := C_Index'Last;
      Initial_Q    : in     A_Matrix := Identity);
 
end Hessenberg;
