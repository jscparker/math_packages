
-- Test QR least squares equation solving real valued square matrices.

with Ada.Numerics.Generic_elementary_functions;
with Givens_QR;
with Test_Matrices;
With Text_IO; use Text_IO;

procedure givens_qr_tst_3 is

   type Real is digits 15;

   subtype Index is Integer range 1..137;

   subtype Row_Index is Index;
   subtype Col_Index is Index;

   Starting_Row : constant Row_Index := Index'First + 0;
   Starting_Col : constant Col_Index := Index'First + 0;
   Final_Row : constant Row_Index := Index'Last - 0;
   Final_Col : constant Col_Index := Index'Last - 0;

   type Matrix is array(Row_Index, Col_Index) of Real;

   type Matrix_inv is array(Col_Index, Row_Index) of Real;
   --  For inverses of A : Matrix; has shape of A_transpose.

   type Matrix_normal is array(Col_Index, Col_Index) of Real;
   --  For A_transpose * A: the Least Squares Normal Equations fit in here.
   --  (A_transpose * A) * x = b are called the normal equations.

   type Matrix_residual is array(Row_Index, Row_Index) of Real;
   --  For A x A_transpose 

 --pragma Convention (Fortran, Matrix); --No! This QR prefers Ada convention.

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;

   package QR is new Givens_QR
     (Real     => Real, 
      R_Index  => Index, 
      C_Index  => Index, 
      A_Matrix => Matrix);
   use QR;

   -- QR exports Row_Vector and Col_Vector

   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;

   package rio is new Float_IO(Real);
   use rio;

   subtype Longer_Real  is Real;     -- general case, and for best speed
 --type Longer_Real  is digits 18;   -- 18 ok on intel, rarely useful

   Min_Real :constant Real := 2.0 ** (Real'Machine_Emin/2 + Real'Machine_Emin/4);

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

 --Desired_Matrix : Matrix_Id  := Upper_Ones;
 --Desired_Matrix : Matrix_Id  := Lower_Ones;
 --Desired_Matrix : Matrix_Id  := All_Ones;           -- singular
 --Desired_Matrix : Matrix_Id  := Zero_Cols_and_Rows; -- singular
 --Desired_Matrix : Matrix_Id  := Easy_Matrix;
 --Desired_Matrix : Matrix_Id  := Symmetric_Banded;
 --Desired_Matrix : Matrix_Id  := Pascal;             -- eigval = x, 1/x   all positive.
 --Desired_Matrix : Matrix_Id  := Forsythe_Symmetric; -- Small diag
 --Desired_Matrix : Matrix_Id  := Forsythe;           -- eigs on circle of rad 2**(-m)
 --Desired_Matrix : Matrix_Id  := Small_Diagonal;
 --Desired_Matrix : Matrix_Id  := Zero_Diagonal_2;
 --Desired_Matrix : Matrix_Id  := Kahan;
 --Desired_Matrix : Matrix_Id  := Non_Diagonalizable; -- one eigval, one eigvec.
 --Desired_Matrix : Matrix_Id  := Peters_0;  -- one small eig., but better than Moler
 --Desired_Matrix : Matrix_Id  := Peters_1;  -- row pivoting likes this
 --Desired_Matrix : Matrix_Id  := Peters_2;  -- max growth with row pivoting
 --Desired_Matrix : Matrix_Id  := Frank; -- eigenvals -> .25; at N=16, 64:1.00...0
 --Desired_Matrix : Matrix_Id  := Ring_Adjacency; -- singular at 4, 8, 12, 16, 24, 32 ..
 --Desired_Matrix : Matrix_Id  := Zero_Diagonal; -- Like Wilkinson, if N odd: 1 zero eig.
 --Desired_Matrix : Matrix_Id  := Wilkinson_minus; -- 1 zero eig N odd, else 1 small eig
 --Desired_Matrix : Matrix_Id  := Moler_0;
 --Desired_Matrix : Matrix_Id  := Moler_1;
 --Desired_Matrix : Matrix_Id  := Random;
 --Desired_Matrix : Matrix_Id  := Vandermonde; -- max size is 256 x 256; else overflows
 --Desired_Matrix : Matrix_Id  := Hilbert;
 --Desired_Matrix : Matrix_Id  := Ding_Dong; -- Eigs clustered near +/- Pi


   A, R, A_lsq : Matrix;
   A_inv_lsq : Matrix_inv;
   A_x_A_inv_minus_I : Matrix_residual;
   A_tr_x_Residual : Matrix_Normal; -- Col_Index x Col_Index
   Q : Q_Matrix;
   Err_in_Least_Squ_A, Frobenius_lsq_soln_err : Real;
   Condition_Number_of_Least_Squ_A, Condition_Number_of_A : Real;
   Sum : Longer_Real;

   Scale   : Col_Vector;
   Permute : Permutation;
   Cutoff_Threshold_1 : constant Real := Two**(-30);
   Singularity_Cutoff_Threshold : Real;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm
     (A            : in     Matrix_normal)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
      return Real
    is
      Max_A_Val : Real := Zero;
      Sum, Scaling, tmp : Real := Zero;
    begin
 
      Max_A_Val := Zero;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Two ** (Real'Machine_Emin + 4);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
      return Real
    is
      Max_A_Val : Real := Zero;
      Sum, Scaling, tmp : Real := Zero;
    begin
 
      Max_A_Val := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Two ** (Real'Machine_Emin + 4);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;


   ------------
   -- Invert --
   ------------
  
   --  Get Inverse of the Matrix:
   
   procedure Get_Inverse 
     (R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      Singularity_Cutoff : in Real;
      A_inverse_lsq  :    out Matrix_inv)     -- shape is A_transpose
   is
      Solution_Row_Vector : Row_Vector; --  X
      Unit_Col_Vector     : Col_Vector := (others => 0.0); --  B
      --  We are going to solve for  X  in  A*X = B

      Final_Row    : Row_Index renames Q.Final_Row;
      Final_Col    : Col_Index renames Q.Final_Col;
      Starting_Row : Row_Index renames Q.Starting_Row;
      Starting_Col : Col_Index renames Q.Starting_Col;
   begin
   
      A_inverse_lsq := (others => (others => Zero));

      --new_line;
      --for i in Starting_Row_Col..Max_Row_Col loop
         --put(R(i,i));
      --end loop;
   
      for i in Starting_Row .. Final_Row loop

         if i > Starting_Row then 
            Unit_Col_Vector(i-1) := 0.0;
         end if;
         Unit_Col_Vector(i) := 1.0;
         --  Make all possible unit Col vectors:  (Final_Row-Starting_Row+1) of them

         QR_Solve 
           (X                  => Solution_Row_Vector,
            B                  => Unit_Col_Vector,
            R                  => R,
            Q                  => Q,
            Row_Scalings       => Scale,
            Col_Permutation    => Permute,
            Singularity_Cutoff => Singularity_Cutoff);
         --  Solve equ. A*Solution_Row_Vector = Unit_Col_Vector (for Solution_Row...).
	 --  Q contains the Starting_Row, Final_Row, Starting_Col, etc.

         for j in Starting_Col .. Final_Col  loop
            A_inverse_lsq (j,i) := Solution_Row_Vector(j);
         end loop;
	 --  All vectors you multiply by A must be Row_Vectors, 
	 --  So the cols of A_inv are Row_Vectors.  

      end loop;
   
   end Get_Inverse;

   procedure Get_Err_in_Least_Squares_A
     (A              : in     Matrix;
      R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      Singularity_Cutoff : in Real;
      A_lsq          :    out Matrix;
      Err_in_Least_Squ_A        : out Real;
      Condition_Number_of_Least_Squ_A : out Real;
      Condition_Number_of_A           : out Real)
   is
      Max_Diag_Val_of_R : constant Real := Abs R(Starting_Col, Starting_Row);
      Min_Allowed_Diag_Val_of_R : constant Real := 
                                      Max_Diag_Val_of_R * Singularity_Cutoff;
      Min_Diag_Val_of_Least_Squ_R : Real;

      Least_Squares_Truncated_Final_Row : Row_Index := Final_Row; -- essential init

      Product_Vector, Col_of_R : Col_Vector := (others => Zero);
      Row : Row_Index;
      Err_Matrix : Matrix := (others => (others => Zero));

      Final_Row    : Row_Index renames Q.Final_Row;
      Final_Col    : Col_Index renames Q.Final_Col;
      Starting_Row : Row_Index renames Q.Starting_Row;
      Starting_Col : Col_Index renames Q.Starting_Col;
   begin

      -- The Columns of R have been permuted; unpermute before comparison of A with Q*R
      -- The Columns of R have been scaled. Before comparison of A with Q*R
      -- Must unscale each col of R by multiplying them with 1/ScalePermute(Col).

      Row := Starting_Row;
      Find_Truncated_Final_Row:
      for Col in Starting_Col .. Final_Col loop
         Min_Diag_Val_of_Least_Squ_R := Abs R(Row, Col);
         if Min_Diag_Val_of_Least_Squ_R < Min_Allowed_Diag_Val_of_R then
	    Least_Squares_Truncated_Final_Row := Row - 1;
            Min_Diag_Val_of_Least_Squ_R := Abs R(Row-1, Col-1);
	    exit Find_Truncated_Final_Row;
	 end if;
         if Row < Final_Row then Row := Row + 1; end if;
      end loop Find_Truncated_Final_Row;

      for Col in Starting_Col .. Final_Col loop

         Col_of_R := (others => Zero);

         for Row in Starting_Row .. Least_Squares_Truncated_Final_Row loop
            Col_of_R(Row) := R(Row, Col);
	 end loop;

         Product_Vector := Q_x_Col_Vector (Q, Col_of_R);

         for Row in Starting_Row .. Final_Row loop
            A_lsq(Row, Permute(Col)) 
	            := Product_Vector(Row) / (Scale(Row) + Min_Real);
            Err_Matrix(Row, Col) :=  A(Row, Permute(Col)) - A_lsq(Row, Permute(Col));
         end loop;
      end loop;

      --  Froebenius norm fractional error = ||Err_Matrix|| / ||A||

      Err_in_Least_Squ_A := 
         Frobenius_Norm (Err_Matrix) / (Frobenius_Norm (A) + Min_Real); 

      Condition_Number_of_Least_Squ_A := 
         Max_Diag_Val_of_R / (Min_Diag_Val_of_Least_Squ_R + Min_Real);

      Condition_Number_of_A := 
         Max_Diag_Val_of_R / (Abs R(Final_Row, Final_Col) + Min_Real);

   end Get_Err_in_Least_Squares_A;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 : String := "") is
     Continue : Character := ' ';
   begin
     New_Line;
     if S0 /= "" then put_line (S0); end if;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     if S8 /= "" then put_line (S8); end if;
     if S9 /= "" then put_line (S9); end if;
     if S10 /= "" then put_line (S10); end if;
     new_line;
     begin
	put ("Type a character to continue: ");
	get_immediate (Continue);
     exception
	when others => null;
     end;
   end Pause;

begin

   Pause( 
    "More on equation solving using the QR decomposition of A. Below we see",
    "if we can find useful solutions for x in A*x = b when the matrix A is",
    "singular or ill-conditioned. If A is singular then there's a vector b such that",
    "no choice of x satisfies A*x - b = 0. If A is not singular but has condition",
    "number 10**40, the desired x exists, but it can have length ~10**40. Below",
    "all equation solving will use the least squares A described in the previous",
    "test. Another way to think of it: A*x = b is solved via R*x = Q'*b. By discarding",
    "Q vectors associated with near-zero diagonal elements of R, we are removing the",
    "component of b that is effectively in the null space of A, and then solving for x.",
    "This doesn't work well unless the decomposition is rank-revealing. That's what we",
    "test below on a test suite of (mostly) ill-conditioned or pathological matrices."
     );


   Pause( 
    "First pass through the matrix test suite: all equation solving will use a least",
    "squares A. Should see: (1) all solutions are good to at least 8 sig. figures,",
    "and (2) no matrix was modified by more than 1 part in 10**9 by the least",
    "squares process. In most cases results are better than that in both categories."
     );


   -- we use Singularity_Cutoff => Two**(-30);  default is Two**(-38)

   for I in 1 .. 2 loop

     Singularity_Cutoff_Threshold := Cutoff_Threshold_1;

     if I > 1 then
       Singularity_Cutoff_Threshold := Zero;
       new_line(3);
       Pause( 
        "Now one last pass through the matrix test suite. In this case there will be",
        "no least squares modification of A.  No Q column vectors are discarded. In",
        "many cases the failure of the equation solving is catastrophic."
        );
     end if;

   for Chosen_Matrix in Matrix_id loop

      Init_Matrix (A, Chosen_Matrix, Index'First, Index'Last);

      R := A; -- cp A to R, then leave A unmolested. R will be overwritten w/ real R
 
      QR_Decompose 
        (A                  => R,  -- input matrix A, but call it R.
         Q                  => Q,
         Row_Scalings       => Scale,
         Col_Permutation    => Permute,
         Final_Row          => Final_Row,
         Final_Col          => Final_Col,
         Starting_Row       => Starting_Row,
         Starting_Col       => Starting_Col);

      Get_Inverse 
        (R, Q, Scale, Permute, Singularity_Cutoff_Threshold, A_inv_lsq);

      Get_Err_in_Least_Squares_A
        (A, R, Q, Scale, Permute, 
	 Singularity_Cutoff_Threshold,
	 A_lsq,
         Err_in_Least_Squ_A, Condition_Number_of_Least_Squ_A, Condition_Number_of_A);

      --  Least Squares solutions of A*x - b = 0 are solutions of the normal
      --  equations    (A_transpose*A) * x - A_transpose*b = 0.
      --  In practice you don't get A*x - b = 0. You find an A*x - b that is
      --  in the null space of A_transpose:  A_transpose * (A*x - b) = 0
      --  
      -- next use QR to solve A*x_i = unit_vec_i over of unit vectors,
      -- The vectors unit_vec_i form the col vecs of I = Identity, and the
      -- vectors x_i form the col vecs of A_Inv. A*A_Inv - Identity /= 0 in
      -- general, but we can test  A_transpose * (A*x - b) = 0, which will
      -- be limited by the condition number of A, and quality of rank revelation.
      --
      -- Get:   A_x_A_inv_minus_I  =   A*A_Inv - Identity
      --
      --  Cols of this mat are the desired residuals. 
      --
      --  A_x_A_inv_minus_I is really   Row_Index x Row_Index:

      for j in Starting_Row .. Final_Row loop
      for i in Starting_Row .. Final_Row loop

         Sum := +0.0;
         for k in Starting_Col .. Final_Col loop
	    Sum := Sum + Longer_Real(A_lsq(i, k)) * Longer_Real(A_inv_lsq(k,j));
         end loop;
	 if i = j then
           A_x_A_inv_minus_I(j, i) := Real (Sum - 1.0);
	 else
           A_x_A_inv_minus_I(j, i) := Real (Sum);
	 end if;
      end loop;
      end loop;

      -- Residual is: A*x - b, (where the x was obtained by solving A*x = b). 
      -- Residual /= 0 unless b is in the space spanned by the col vectors of A.
      -- Residual won't be 0 in general unless A is full rank (non-singular),
      -- (so that its cols span all space).  This may be rare!
      -- Remember, our least squares A_lsq, obtained by discarding Q vecs is 
      -- is not actually full rank; A_lsq_inverse * A_lsq /= I. However if you
      -- restrict the calculation to the space spanned by the col vectors of A
      -- (multiply on the left by A_tr), can use A_tr*(A_lsq_inverse * A_lsq - I)
      -- as test of success of equation solving.
      -- So we get
      --
      --     A_tr_x_Residual = A_transpose * (A * A_inv - I)
      --
      --  result has shape: Col_Index x Col_Index

      for j in Starting_Col .. Final_Col loop
      for i in Starting_Col .. Final_Col loop

         Sum := +0.0;
         for k in Starting_Col .. Final_Col loop
	    Sum := Sum + Longer_Real(A_lsq(k, i)) * Longer_Real(A_x_A_inv_minus_I(k,j));
         end loop;
         A_tr_x_Residual(j, i) := Real (Sum);
      end loop;
      end loop;

      --   fractional error:
      --       = |A_tr*Residual| / |A_tr| =  |A_tr*(Ax - b)| / |A_tr|

      --   Frobenius_lsq_soln_err  
      --    := Max_Col_Length (A_tr_x_Residual) / (Max_Row_Length (A_lsq) + Min_Real);
      --
      Frobenius_lsq_soln_err  
          := Frobenius_Norm (A_tr_x_Residual) / (Frobenius_Norm (A_lsq) + Min_Real);

      new_line;
      put("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));  put(":");
      new_line;
      put(" Approx condition number    (lower bound) of original  A =");
      put(Condition_Number_of_A);
      new_line;
      put(" Approx condition number    (lower bound) of least squ A =");
      put(Condition_Number_of_Least_Squ_A);
      new_line;
      put(" Approx. of max|A'*(A*x-b)| / |A'*b| over unit vectors b ="); 
      put(Frobenius_lsq_soln_err);
      new_line(1);
      -- actually, best to use a different norm to get max|A'*(A*x-b)| / |A'*b| over b's,
      -- (max col length probably) but some other time maybe.

   end loop;
   end loop;

end;
