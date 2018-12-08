
---------------------------------------------------------------------------
-- package QR_Symmetric_Eigen, QR based eigen-decomposition
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

-- package QR_Symmetric_Eigen
--
-- Eigen-decomposition of symmetric real-valued matrices.
-- 
-- The matrix is tri-diagonalized with Givens rotations, then diagonalized
-- using QR iterations.
-- 
-- procedure Eigen_Decompose 
--
-- Works on arbitrary diagonal blocks of input matrix. For other blocks just
-- copy the matrix to desired position; copy overhead is small compared
-- to the O(N^3) running time of the decomposition.

generic

   type Real is digits <>;
   type Index is range <>;
   type Matrix is array (Index, Index) of Real;

package QR_Symmetric_Eigen is

   type Col_Vector is array(Index) of Real;

   -- procedure Eigen_Decompose 
   --
   -- The routine returns eigenvectors and eigenvalues of arbitrary diagonal
   -- blocks of any real-valued square symmetric matrix.
   --
   -- The orthonormal (unordered) eigenvectors are the Columns of Q.
   -- The orthonormal (unordered) eigenvectors are returned as the Rows of Q'=Q_tr.
   -- Eigenvals (returned in array Eigenvals) are ordered the same as Eigvecs in Q.
   -- So A = Q * E * Q'. The diagonal elements of diagonal matrix E are the
   -- eigvals.  The routine performs the eigen-decomposition on arbitrary square
   -- diagonal blocks of matrix A. 
   -- It is assumed that the blocks are symmetric matrices.
   -- The upper left corner of the square matrix is (Start_Col, Start_Col).
   -- The lower rgt  corner of the square matrix is (Final_Col, Final_Col).
   -- Matrix A doesn't need to be positive definite, or semi-definite.
   -- If Eigenvectors_Desired = False, then Q is not calculated.
   --
   -- Input matrix A is destroyed. Save a copy of A if you need it.
   --
   procedure Eigen_Decompose
     (A                      : in out Matrix;     -- destroyed
      Q                      :    out Matrix;     -- columns of Q are the eigvecs
      Eigenvals              :    out Col_Vector;
      Start_Col              : in     Index        := Index'First;
      Final_Col              : in     Index        := Index'Last;
      Eigenvectors_Desired   : in     Boolean      := False);

   procedure Sort_Eigs
     (Eigenvals         : in out Col_Vector;
      Q                 : in out Matrix;          -- columns of Q are the eigvecs
      Start_Col         : in     Index   := Index'First;
      Final_Col         : in     Index   := Index'Last;
      Sort_Eigvecs_Also : in     Boolean := False);

   function Norm
     (Q            : in Col_Vector;
      Starting_Col : in Index := Index'First;
      Final_Col    : in Index := Index'Last)
      return Real;

end QR_Symmetric_Eigen;

