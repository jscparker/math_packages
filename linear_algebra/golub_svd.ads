
---------------------------------------------------------------------------
-- package Golub_SVD, Singular Value Decomposition
-- Copyright (C) 2018 Jonathan S. Parker.
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

--  package Golub_SVD
--
--  Golub-Reinsch Singular Value Decomposition.
--
--  A rectangular real-valued matrix A is bidiagonalized using
--  Householder reflections, followed by QR iterations to complete
--  the diagonalization. A is decomposed into A = U*W*V' where the
--  diagonal of matrix W contains the singular vals (it's really just
--  a 1-dim array here), and where U and V are orthogonal: U'*U = I
--  and V'*V = I.  If A is an m x n matrix (m rows, n cols) then
--  then V is n x n and U is m x m. W is conceptually m x n.
--
--  The eigenvalues of A'A are the squares of the singular vals.
--  If A = U*W*V', U'*U = I and V'*V = I then A'*A = V*W'*U'*U*W*V',
--  so that A'*A * V = V * W'*W. Matrix W is a diagonal matrix with
--  the singular values down the diagonal, so W'*W has the squares
--  of the singular values down its diagonal. In other words the col
--  vectors of V are the eigenvectors of A'A, with eigenvalues given
--  by the diagonal of W'*W.  
--
--  Similarly, A*A' = UWV'VW'U' so that A*A'* U = U * W*W' (using
--  V'V = I). In other words the columns of U are the eigenvectors
--  of A*A'.
--
--  The code is based on the LINPACK Fortran 90 SVD, with several
--  changes to improve reliability and numerical precision.
--
--  All of the Givens rotations were enhanced to improve numerical 
--  precision. Givens rotations are sluggish, but numerical precision 
--  is favored over speed here. The enhanced precision rotations also
--  improved the quality of the orthogonal matrices U and V. Loss in
--  speed is perceptible, but not huge (< 15%).
--  
--  Notes on use
--  ------------
--  
--  You instantiate the generic with Row_Index and Col_Index satisfying
--
--     Row_Index'First = Col_Index'First
--
--  The routine forces Starting_Row = Starting_Col. In other words
--  you input Starting_Col, and Starting_Row is chosen for you.
--  
--  SVD works on square matrices, or matrices with more rows than
--  columns. If the input matrix A has more columns than rows, then
--  let B = transpose(A), and SVD matrix B to get B = U*W*V'. Take
--  the transpose of U*W*V' to get the SVD of matrix A. In other
--  words A = V*W*U'.
--  
--  The procedure works on rectangular sub-blocks of the input matrix A,
--  but the sub-block must have a corner on the diagonal: SVD_Decompose 
--  assumes that the id of the Starting Column equals the id of the
--  Starting Row. So the procedure is convenient for doing diagonal blocks
--  of a matrix, but if you want to SVD arbitrary sub-matrices of the
--  matrix, then you must copy the submatrix to one satisfying the above 
--  constraint. (That's the usual way in these routines anyway; the
--  overhead is negligible in comparison to the SVD.)
--
generic

  type Real is digits <>;

  type Row_Index is range <>;
  type Col_Index is range <>;
  --  You can use different types for Row_Index and Col_Index,
  --  but it doesn't make much sense.  Compiler might find
  --  it easier to optimize if they are both subtypes of Integer.

  type A_Matrix is array (Row_Index, Col_Index) of Real;
  --  Place pragma Convention (Fortran, A_Matrix) just after 
  --  declaration of type A_Matrix. It matters!

package Golub_SVD is

  subtype U_Col_Index is Row_Index; 
  --  If you choose Row_Index here then U is square. This is the stnd choice.

  type U_Matrix is array (Row_Index, U_Col_Index) of Real;
  pragma Convention (Fortran, U_Matrix); 

  type V_Matrix is array (Col_Index, Col_Index) of Real;
  pragma Convention (Fortran, V_Matrix); 

  type Singular_Vector is array (Col_Index) of Real;

  procedure SVD_Decompose (
     A                : in out A_Matrix;         --  m x n  (A is destroyed)
     U                :    out U_matrix;         --  m x m
     V                :    out V_matrix;         --  n x n
     S                :    out Singular_Vector;
     Id_of_1st_S_Val  :    out Col_Index;        --  see 1st remark below
     Starting_Col     : in     Col_Index;        -- Starting_Row is set to Starting_Col
     Final_Col        : in     Col_Index;
     Final_Row        : in     Row_Index;
     Matrix_U_Desired : in     Boolean   := True;
     Matrix_V_Desired : in     Boolean   := True);

     --  If all singular values were judged to be successfully calculated, then
     --     Id_of_1st_S_Val  =  Starting_Col
     --
     --  QR convergence is judged successful for:
     --     S(Id_of_1st_S_Val), S(Id_of_1st_S_Val+1), S(Id_of_1st_S_Val+2), etc.
     --
     --  SVD_Decompose requires: at least 2 row and 2 cols.
     --
     --  SVD_Decompose requires: Integer(Col_Index'First) = Integer(Row_Index'First)
     --    (only because Starting_Row is set to Starting_Col; get unexpected
     --     failures if this is not enforced from the beginning.)
     --
     --  If that's not OK, then copy the matrix to one that satisfies
     --  this constraint. Overhead from this O(N^2) copy is usually negligible.

  Max_Allowed_No_of_Iterations : constant := 2**24;
  --  The QR algorithm iterates for the Singular Vals, but convergence is not
  --  guaranteed. On the other hand, some matrices (eg Vandermode, Kahan, Pascal)
  --  use very large numbers of iterations here, but still finish very fast.
  --  Large numbers of QR iterations seem to occur when the matrix has bad 
  --  row scaling.  The above standard setting (2**24) is ok. (16000x16000
  --  Vandermode matrices use about 2**13 iterations, but take same time to
  --  SVD as ordinary matrices.)

  SVD_Convergence_Error : exception;  
  --  Raised when nothing at all seems to converge: in other words none of the
  --  calculated singular vals are expected to be reliable.

end Golub_SVD;


