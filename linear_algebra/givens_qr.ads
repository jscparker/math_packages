
--------------------------------------------------------------------------
-- package Givens_QR, QR decomposition using Givens rotations
-- Copyright (C) 2008-2018 Jonathan S. Parker
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

--  package Givens_QR
--  
--  QR decomposition and linear equation solving for rectangular real
--  valued matrices. The procedure uses Givens rotations with column 
--  pivoting and row scaling to factor matrix A into a product of an
--  orthogonal matrix Q and an upper triangular matrix R:  A*S*P = Q*R. 
--  Matrix S is a diagonal matrix that scales the columns of A. Matrix
--  P contains the column permutations performed during col pivoting.
--  
--  A matrix with m rows, n columns (m x n) and elements of generic type
--  Real is input as "A" and returned in Q R form.  Matrix A is overwritten
--  with R and should be saved prior to decomposition if it needed. The
--  QR decomposition can be performed on arbitrary sub-blocks of matrix A.
--  Must have No_of_Rows >= No_Of_Cols: m >= n.  (If m < n you take the QR
--  decomposition of A_transpose to get an LQ decomposition of A; see below.)
--  
--  The decomposition uses Givens rotations rather than the more common
--  Householder reflections. In the present package, the Givens arithmetic
--  is enhanced to provide a bit more numerical accuracy.
--
--  Column pivoting is slow, but it's a necessary condition for calculating
--  reliable least squares solutions of x in A*x = b, and makes QR applicable
--  in a lot of cases that would otherwise require an even slower singular
--  value decomposition (SVD).
--  
--  Procedures QR_Decompose and QR_Solve return solutions of the simultaneous 
--  linear equations A*x = b, where A is an arbitrary real valued matrix.
--  The vector b is input into procedure QR_Solve, and the solution is 
--  returned as x. QR can be used to generate least squares 
--  solutions: vectors x that minimize the sum of squares: |A*x - b|^2.
--  (The system of linear equations A*x = b may have no exact solution x, or
--  it may have an infinite number of solutions.)
--  
--  The QR assumes m >= n. (Number of Rows is greater than or equal to 
--  number of Cols.)  If m < n (an underdetermined system of equations)
--  then you QR decompose the transpose of A.  The user of the routine
--  does the transpose, not the routine.  (See Notes on Use below.)
--  
--  Notes on Use
--  
--  Notice we require No_of_Rows >= No_Of_Cols: m >= n.
--
--   A*x = b is
--  
--   a  a  a  a        x       b
--   a  a  a  a        x       b
--   a  a  a  a    *   x   =   b
--   a  a  a  a        x       b
--   a  a  a  a                b
--   a  a  a  a                b
--   a  a  a  a                b
--   a  a  a  a                b
--  
--  Suppose instead that m < n.  (Number of rows is less than number of 
--  columns.)  Then  A x = b is an underdetermined system of
--  equations.  In this case there may exist an infinite
--  number of solutions for x. Usually the preferred solution x is the
--  x of minimum length. To solve for the x of minimum length, decompose 
--  the transpose of A. Let A' = transpose(A) and obtain the usual Q*R = A'.  
--  Using Q and R in this form yields solutions of x in A x = b that have the
--  desired property. Here the equation A x = b becomes R' (Q' x) = b, 
--  or R' z = b, where R' is m x m, and z and b are vectors of length m.

generic

   type Real is digits <>;

   type R_Index is range <>;  --  Row_id : R_Index; 
   type C_Index is range <>;  --  Col_id : C_Index;
   --  This means:
   --  Col vectors are indexed by R_Index, 
   --  Row vectors are indexed by C_Index. 
   --  The matrix is M x N:  M rows by N cols.
   --  M rows means that:    M = R_Index'Last - R_Index'First + 1
   --  N cols means that:    N = C_Index'Last - C_Index'First + 1

   type A_Matrix is array(R_Index, C_Index) of Real;
   --  A is the matrix to be QR decomposed.
   --  The Q matrix (a sequence of Givens rotations) is stored in a 
   --  matrix with the same shape as A.


package Givens_QR is

   type Row_Vector is array(C_Index) of Real; 
   --  Row vector of A_Matrix and R_Matrix. (R is essentially a square matrix)

   type Col_Vector is array(R_Index) of Real; 
   --  Col vector of A_Matrix

   type Permutation is array(C_Index) of C_Index;
   --  The columns are rearranged as a result of pivoting. The swaps are
   --  stored in Col_Permutation : Permutation.  So if Col 1 is swapped with
   --  Col 6, then Col_Permutation(1) = 6, and Col_Permutation(6) = 1.

   type Rotation is record
     Cosine : Real;
     Sine   : Real;
   end record;

   Identity : constant Rotation := (+0.0, +0.0);

   type Rot_Set  is array (R_Index, C_index) of Rotation;
   type Rot_Code is array (R_Index, C_index) of Boolean;
   pragma Pack (Rot_Code);

   type Q_matrix is record
      Rot              : Rot_Set;
      P_bigger_than_L  : Rot_Code;
      Final_Row    :  R_Index;
      Final_Col    :  C_Index;
      Starting_Row :  R_Index;
      Starting_Col :  C_Index;
   end record;

   --  Notice   Q_matrix.Rot   has the same shape as A: M x N.
   --  The true Q is M x M, not M x N like A, but the above matrix holds
   --  the rotation matrices that can be used to create Q.


   procedure QR_Decompose 
     (A                  : in out A_Matrix;     --  overwritten with R
      Q                  :    out Q_Matrix;
      Row_Scalings       :    out Col_Vector;
      Col_Permutation    :    out Permutation;
      Final_Row          : in     R_Index := R_Index'Last;
      Final_Col          : in     C_Index := C_Index'Last;
      Starting_Row       : in     R_Index := R_Index'First;
      Starting_Col       : in     C_Index := C_Index'First);
   
   --  procedure QR_Decompose factors A into  A*S*P = Q*R, where S is diagonal,
   --  and P permutes the order of the columns of a matrix.
   --  
   --  The input matrix A is overwritten with R.
   --  
   --  Must have No_of_Rows >= No_Of_Cols, and at least 2 rows and at least 2 cols.
   --  otherwise Argument_Error is raised:
   --  
   --  if No_Of_Cols > No_of_Rows then   raise Ada.Numerics.Argument_Error;  end if;
   --  if No_Of_Cols < 2.0        then   raise Ada.Numerics.Argument_Error;  end if;


   Default_Singularity_Cutoff  : constant Real := Real'Epsilon * (+2.0)**12;
   --  Used by QR_Solve to handle singular matrices, and ill-conditioned matrices:
   --  
   --  If a diagonal element of R falls below this fraction of the max R element,
   --  then QR_Solve tosses out the corresponding Q column (to get a matrix that
   --  is very near the original matrix, but much better conditioned). The range
   --  Real'Epsilon * 2**10 to Real'Epsilon * 2**24 seems ok. Ill-conditioned
   --  matrices can behave either worse or better as you vary over this range. 
   --  2.0**12 is a good compromise. It means Default_Singularity_Cutoff ~ 2**(-38).

   procedure QR_Solve 
     (X                  :    out Row_Vector;    --  Solution vector
      B                  : in     Col_Vector;
      R                  : in     A_Matrix; --  Original A holds the R after call to QR
      Q                  : in     Q_Matrix; 
      Row_Scalings       : in     Col_Vector;
      Col_Permutation    : in     Permutation;
      Singularity_Cutoff : in     Real := Default_Singularity_Cutoff);


   --  Get Product = Q'*B

   function Q_transpose_x_Col_Vector
     (Q : in Q_Matrix;
      B : in Col_Vector)
      return Col_Vector;


   --  Get Product = Q*B   

   function Q_x_Col_Vector
     (Q : in Q_Matrix;
      B : in Col_Vector)
      return Col_Vector;

   --  type A_Matrix is a list of Col_Vector's:

   procedure Q_x_Matrix
     (Q : in Q_Matrix;
      A : in out A_Matrix);


  --  type V_Index is range <>;

  --  Get Product = Q*V

   subtype V_Index is R_Index;  
   --  This makes V is square, so can be used for making Q in 
   --  explicit matrix form.

   type V_Matrix is array(R_Index, V_Index) of Real;
   --  V might for example be the identity matrix so can get Q in 
   --  matrix form: Q*I = Q, using    procedure Q_x_V_Matrix.
 
   procedure Q_x_V_Matrix
    (Q : in     Q_Matrix;
     V : in out V_Matrix);


   procedure Q_transpose_x_V_Matrix
    (Q : in     Q_Matrix;
     V : in out V_Matrix);

end Givens_QR;
