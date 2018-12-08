
-- Test Jacobi Eigendecomposition of real valued square matrices.

with Ada.Numerics.Generic_elementary_functions;
with Jacobi_Eigen;
with Test_Matrices;
with Text_IO; use Text_IO;

procedure jacobi_eigen_tst_1 is

   type Real is digits 15;

 --Matrix_Size : constant := 2277;
 --Matrix_Size : constant := 597;
   Matrix_Size : constant := 137;

   Index_First : constant := 1;

   -- Sometimes it's faster if you use a storage matrix that's a little
   -- larger than the original matrix. (Depends on the hardware, the problem,
   -- and the size of the matrix.) A rule of thumb: Width of storage matrix
   -- should be an even number, greater than matrix size, and never 2**n.
   --
   -- Most important thing is to avoid 2**n matrix width.
   --
   -- Padding usually doesn't help a lot, but just in case let's add 1 or 2:

   Padding : constant := 2 - (Matrix_Size - Index_First + 1) mod 2;

   subtype Index is Integer range 1 .. Matrix_Size + Padding;

   -- the test matrix is square-shaped matrix on:  Index x Index.
   -- eg Hilbert's matrix is a square matrix with unique elements on the range
   -- Index'First .. Index'Last.  However, you  have the option or using any 
   -- diagonal sub-block of the matrix defined by Index x Index

   subtype Col_Index is Index;

   Starting_Col : constant Index := Index'First + 0;
   Final_Col    : constant Index := Index'Last  - Padding;

   -- Can't change:

   Starting_Row : constant Index := Starting_Col;
   Final_Row    : constant Index := Final_Col;

   type Matrix is array(Index, Index) of Real;

   --pragma Convention (Fortran, Matrix); --No! prefers Ada convention.

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;

   package Eig is new Jacobi_Eigen
     (Real   => Real, 
      Index  => Index, 
      Matrix => Matrix);
   use Eig;

   -- Eig exports Col_Vector

   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;

   package rio is new Float_IO(Real);
   use rio;

 --subtype Real_Extended is Real;       -- general case, works fine
   type Real_Extended is digits 18;   -- 18 ok on intel

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

   A, A_true, Q_tr : Matrix;
   Eigenvals : Col_Vector;
   Frobenius_QtrQ_Err, Frobenius_QQtr_Err, Frobenius_QEQ_Err : Real;

   No_of_Sweeps_Done, No_of_Rotations : Natural;
   
   N : constant Real := Real (Starting_Col) - Real (Final_Col) + 1.0;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A : in Matrix)
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

   ------------------------------------
   -- Get_Err_in_Reassembled_Q_and_A --
   ------------------------------------
  
   --  check that A = V*E*Q_tr
   --  E is diagonal with the eigs along the diag.
   --  V is orthogonal with the eig vecs as columns.

   
   procedure Get_Err_in_Reassembled_Q_and_A 
     (A                  : in     Matrix;  -- true original A
      Q_tr               : in     Matrix;
      E                  : in     Col_Vector;
      Final_Col          : in     Col_Index;
      Starting_Col       : in     Col_Index;
      Frobenius_QtrQ_Err :    out Real;
      Frobenius_QQtr_Err :    out Real;
      Frobenius_QEQ_Err  :    out Real)
   is
      Err, S          : Real;
      Min_Real : constant Real := +2.0 **(Real'Machine_Emin + 4);
    
      Sum : Real_Extended;

      Identity    : Matrix := (others => (others => Zero));
      Product_QQ  : Matrix := (others => (others => Zero));
      Product_QEQ : Matrix := (others => (others => Zero));

      subtype Index_Subrange is Index range Starting_Col .. Final_Col;

   begin
  
      for r in Index_Subrange loop
        Identity(r, r) := 1.0;
      end loop;

      -- Find error in I - Q*Q' etc.
      -- Notation: Q' == Q_tr == transpose of Q.

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := +0.0;
        for j in Index_Subrange loop
           Sum := Sum + Real_Extended (Q_tr(j, Row)) * Real_Extended (Q_tr(j, Col)); 
        end loop;
        Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- Get Frobenius norm of:  Product_QQ - I:

      S := Zero;
      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
         Err := Identity(Row, Col) - Product_QQ(Row, Col);
         S := S + Err * Err;
      end loop;
      end loop;

      -- Get fractional Frobenius :  Frobenius (Product_QQ - I) / Frobenius (I):

      Frobenius_QQtr_Err := Sqrt(S) / Sqrt (-Real(Starting_Col)+Real(Final_Col)+1.0);
  

      -- Find error in I - Q'*Q.
      -- reuse array Product_QQ:

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := +0.0;
        for j in Index_Subrange loop
           Sum := Sum + Real_Extended (Q_tr(Row, j)) * Real_Extended (Q_tr(Col, j)); 
        end loop;
        Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- Get Frobenius norm of:  Product_QQ - I:

      S := Zero;
      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
         Err := Identity(Row, Col) - Product_QQ(Row, Col);
         S := S + Err * Err;
      end loop;
      end loop;

      -- Get fractional Frobenius :  Frobenius (Product_QQ - I) / Frobenius (I):

      Frobenius_QtrQ_Err := Sqrt(S) / Sqrt (-Real(Starting_Col)+Real(Final_Col)+1.0);
   
      --  check that A = Q*E*Q_tr
      --  E is diagonal with the eigs along the diag.
      --  Q is orthogonal with the eig vecs as columns.

      -- explicitly calculate Q*E*Q_tr:

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := +0.0;
        for j in Index_Subrange loop
           Sum := Sum +
              Real_Extended (Q_tr(j, Row)) *           --  Q(Row, j)
              Real_Extended (E(j))         *           --  j-th eig is const along Q col
              Real_Extended (Q_tr(j, Col)); 
        end loop;
        Product_QEQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- resuse array Product_QEQ to get Error Matrix := Product_QEQ - A:

      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
         Product_QEQ(Row, Col) := A(Row, Col) - Product_QEQ(Row, Col);
      end loop;
      end loop;

      Frobenius_QEQ_Err := Frobenius_Norm (Product_QEQ) / 
                                   (Frobenius_Norm (A) + Min_Real);
   
   end Get_Err_in_Reassembled_Q_and_A;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 : string := "") is
     Continue : Character := ' ';
   begin
     new_line;
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
     new_line;
     begin
        put ("Type a character to continue: ");
        get_immediate (Continue);
     exception
        when others => null;
     end;
   end pause;

begin

   Pause( 
     "Test 1: Jacobi Eigendecomposition of matrix A.",
     " ",
     "The Jacobi Eigendecomposition of A is successful if the identities Q*Q' = I,",
     "Q'*Q = I and Q*E*Q' = A are satisfied. Here Q' denotes the transpose of Q, and",
     "E is any diagonal matrix. (E will hold the Eigvals.) If 15 digit Reals are used,",
     "then we expect the error in the calculation of A = Q*E*Q' to be (hopefully)",
     "a few parts per 10**15. In other words ||Q*E*Q' - A|| / ||A|| should be a few",
     "multiples of 10**(-15). Here ||*|| denotes the Frobenius Norm. Other matrix",
     "norms give slightly different answers, so its an order of magnitude estimate."
     );

   --  Get A = Q*E*Q_tr, or E = Q_tr*A*Q.  
   --  E is diagonal with the eigs along the diag.
   --  V is orthogonal with the eig vecs as columns.

   for Chosen_Matrix in Matrix_Id loop

      Init_Matrix (A, Chosen_Matrix, Starting_Col, Final_Col);

      --  Usually A is not symmetric.  Eigen_Decompose doesn't care about
      --  that. It uses the upper triangle of A, and pretends that A is
      --  symmetric. But for subsequent analysis,  we want it symmetrized:

      for c in Starting_Col .. Final_Col  loop
      for r in Starting_Row .. Final_Row  loop
         A_true(r, c) := 0.5 * (A(c, r) + A(r, c));
      end loop;
      end loop;

      A := A_true; -- A will be destroyed, so save A_true
 
      Eigen_Decompose 
        (A                      => A,   -- A is destroyed
         Q_tr                   => Q_tr,
         Eigenvals              => Eigenvals,
         No_of_Sweeps_Performed => No_of_Sweeps_Done,
         Total_No_of_Rotations  => No_of_Rotations,
         Final_Col              => Final_Col,
         Start_Col              => Starting_Col,
         Eigenvectors_Desired   => True);

      Sort_Eigs
        (Eigenvals          => Eigenvals,
         Q_tr               => Q_tr,
         Start_Col          => Starting_Col,
         Final_Col          => Final_Col,
         Sort_Eigvecs_Also  => True);

    --put(Q_tr(1,1));

      if False then
      for I in Starting_Col .. Final_Col loop
         put(Eigenvals(I)); 
      end loop;
      new_line(2);
      end if;


      new_line;
      put("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));

      Get_Err_in_Reassembled_Q_and_A 
        (A                 => A_True, 
         Q_tr              => Q_tr, 
         E                 => Eigenvals,
         Final_Col         => Final_Col, 
         Starting_Col      => Starting_Col, 
         Frobenius_QtrQ_Err => Frobenius_QtrQ_Err, 
         Frobenius_QQtr_Err => Frobenius_QQtr_Err, 
         Frobenius_QEQ_Err  => Frobenius_QEQ_Err);

      --  Froebenius norm fractional error:
      --                Max_Error_F  = ||Err_Matrix|| / ||A||

      new_line;
      put(" No of sweeps performed, and Total_No_of_Rotations / (N*(N-1)/2) ="); 
      new_line;
      put(Real (No_of_Sweeps_Done)); 
      put(Real (No_of_Rotations) / (N*(N-1.0)/2.0));
      new_line;
      put(" Err in I-Q*Q'   (Q = orthogonal) is  ||I-Q*Q'|| / ||I|| ="); 
      put(Frobenius_QQtr_Err);
      new_line;
      put(" Err in I-Q'*Q   (Q = orthogonal) is  ||I-Q'*Q|| / ||I|| ="); 
      put(Frobenius_QtrQ_Err);
      new_line;
      put(" Err in A-Q*E*Q' (E = eigenvals) is ||A-Q*E*Q'|| / ||A|| ="); 
      put(Frobenius_QEQ_Err);
      new_line;

      --Pause; new_line;

   end loop;

end jacobi_eigen_tst_1;
