
-- Test Jacobi Eigendecomposition of real valued square matrices.
-- Uses Intel's 18 digit Real to estimate eigval err in 15 digit real,
-- so not very portable.

with Jacobi_Eigen;
with Test_Matrices;
with Text_IO; use Text_IO;

procedure jacobi_eigen_tst_2 is

   type Index is range 1 .. 138;

   -- the test matrix is square-shaped matrix on:  Index x Index.
   -- eg Hilbert's matrix is a square matrix with unique elements on the range
   -- Index'First .. Index'Last.  However, you  have the option or using any 
   -- diagonal sub-block of the matrix defined by Index x Index

   Starting_Col : constant Index := Index'First + 0;
   Final_Col    : constant Index := Index'Last  - 0;

   -- Can't change:

   Final_Row    : constant Index := Final_Col;


   type Real is digits 18;

   type Matrix is array(Index, Index) of Real;

   --pragma Convention (Fortran, Matrix); --No! prefers Ada convention.

   package Eig is new Jacobi_Eigen
     (Real   => Real, 
      Index  => Index, 
      Matrix => Matrix);
   use Eig;

   -- Eig exports Col_Vector

   package rio is new Float_IO(Real);
   use rio;

   type Real_15 is digits 15;

   type Matrix_15 is array(Index, Index) of Real_15;

   A : Matrix_15;

   --pragma Convention (Fortran, Matrix); --No! prefers Ada convention.

   package Eig_15 is new Jacobi_Eigen
     (Real   => Real_15, 
      Index  => Index, 
      Matrix => Matrix_15);

   -- Eig exports Eig_15.Col_Vector

   package Make_Square_Matrix is new Test_Matrices (Real_15, Index, Matrix_15);
   use Make_Square_Matrix; 

   Eig_vals_15 : Col_Vector := (others => 0.0);

  procedure d15 
    (A      : in  Matrix_15;
     Eigval : out Col_Vector) is

   Q_tr, A_tmp : Matrix_15;
   Eigenvalues : Eig_15.Col_Vector := (others => 0.0);

   No_of_Sweeps_Done, No_of_Rotations : Natural;
   
  begin
      A_tmp := A;

      Eig_15.Eigen_Decompose 
        (A                      => A_tmp,   -- A_tmp is destroyed
         Q_tr                   => Q_tr,
         Eigenvals              => Eigenvalues,
         No_of_Sweeps_Performed => No_of_Sweeps_Done,
         Total_No_of_Rotations  => No_of_Rotations,
         Final_Col              => Final_Col,
         Start_Col              => Starting_Col,
         Eigenvectors_Desired   => True);

      Eig_15.Sort_Eigs
        (Eigenvals          => Eigenvalues,
         Q_tr               => Q_tr,
         Start_Col          => Starting_Col,
         Final_Col          => Final_Col,
         Sort_Eigvecs_Also  => True);

        for i in Starting_Col .. Final_Col  loop
            Eigval(i) := Real (Eigenvalues(i));
        end loop;

  end d15;


   Zero : constant Real := +0.0;

   B, Q_tr : Matrix := (others => (others => Zero));
   Eigenvals : Col_Vector;
   Err, Ave_Err, Max_Err, Max_Val : Real;
   Min_Allowed_Real : constant Real := +2.0 **(Real'Machine_Emin + 32);

   No_of_Sweeps_Done, No_of_Rotations : Natural;
   
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

   Pause (
      "The test uses Intel's 18 digit Reals to estimate eigenval error in the",
      "15 digit calculation, so only works on Intel/AMD hardware."
         );

   for Pass_id in 1 .. 2 loop
   for Chosen_Matrix in Matrix_Id loop

      Init_Matrix (A, Chosen_Matrix, Starting_Col, Final_Col);

      --  Usually A is not symmetric.  Eigen_Decompose doesn't care about
      --  that. It uses the upper triangle of A, and pretends that A is
      --  symmetric. But for subsequent analysis,  we symmetrize:

      if Pass_id = 1 then -- use lower triangle of A to make a fully symmetric A:

        for Col in Starting_Col .. Final_Col  loop
        for Row in Col .. Final_Row  loop
           A(Col, Row) := A(Row, Col);  -- write lower triangle of A to upper triangle
        end loop;
        end loop;

      else                -- use upper triangle of A to make a fully symmetric A:

        for Col in Starting_Col .. Final_Col  loop
        for Row in Starting_Col .. Col  loop
           A(Col, Row) := A(Row, Col);  -- write lower triangle of A to upper triangle
        end loop;
        end loop;

      end if;

      d15 (A, Eig_vals_15);

      for i in Starting_Col .. Final_Col  loop
      for j in Starting_Col .. Final_Row  loop
         B(i,j) := Real (A(i,j)); 
      end loop;
      end loop;

      Eigen_Decompose 
        (A                      => B, 
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

      Max_Err := Zero;
      Ave_Err := Zero;
      for I in Starting_Col .. Final_Col loop
         Err := Abs (Eigenvals(I) - Eig_vals_15(I)); 
         if  Err > Max_Err  then  Max_Err := Err;  end if;
         Ave_Err := Ave_Err + Err;
      end loop;

      Max_Val := Zero;
      for I in Starting_Col .. Final_Col loop
         if Abs (Eigenvals(I)) > Max_Val then Max_Val := Abs (Eigenvals(I));  end if;
      end loop;

      new_line;
      put("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));
      new_line;
      put(" Max err:");
      put(Max_Err / (Max_Val + Min_Allowed_Real));
      put("  Ave err:");
      put(Ave_Err / ((Real(Final_Col)-Real(Starting_Col)+1.0)*(Max_Val+Min_Allowed_Real)));

   end loop;
   end loop;

end;
