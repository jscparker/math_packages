
-- Test Eigendecomposition of Real valued square matrices.

with Test_Matrices;
with QR_Symmetric_Eigen;
with Text_IO; use Text_IO;

procedure qr_symmetric_eigen_tst_1 is

   type Real is digits 15;
 --type Real is digits 18;

 --Matrix_Size : constant := 7000;
 --Matrix_Size : constant := 1303;
   Matrix_Size : constant := 137;
 --Matrix_Size : constant := 57;
 
   Padding : constant := 1;

   subtype Index is Integer range 1 .. Matrix_Size + Padding;

   -- the test matrix is square-shaped matrix on:  Index x Index.

   -- Can't change:

   Starting_Col : constant Index := Index'First;
   Final_Col    : constant Index := Index'Last - Padding;

   subtype Submatrix_Range is Index range Starting_Col .. Final_Col;

   type Matrix is array(Index, Index) of Real;

   pragma Convention (fortran, Matrix);

   package QR_Method is new QR_Symmetric_Eigen (Real, Index, Matrix);
   use QR_Method;

   A        : Matrix;
   Eigvals  : Col_Vector;

   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;

   package rio is new Float_IO(Real);
   use rio;

 --subtype Real_e is Real;     -- general case, works fine
   type Real_e is digits 18;   -- 18 ok on intel

   type Col_Vector_e is array (Index) of Real_e;

   Min_Exp          : constant Integer := Real'Machine_Emin - Real'Machine_Emin/4;
   Min_Allowed_Real : constant Real    := 2.0**Min_Exp;  

   -------------
   -- Product --
   -------------
  
   function Product
     (A            : in Matrix;
      Vector       : in Col_Vector_e;
      Starting_Col : in Index;
      Final_Col    : in Index)
      return Col_Vector_e is
      Result : Col_Vector_e := (others => 0.0);
      Sum : Real_e;
   begin

      for Row in Starting_Col .. Final_Col loop
         Sum := +0.0;
         for c in Starting_Col .. Final_Col loop
            Sum := Sum + Real_e (A(Row, c)) * Vector(c);
         end loop;
         Result (Row) := Sum;
      end loop;

      return Result;

   end Product;

   ----------------------------
   -- Error_in_Decomposition --
   ----------------------------

   -- rough estimate, relative err.
   -- most of this only works for A'Range, not sub range.

   procedure Error_in_Decomposition
     (A            : in     Matrix;
      W_r          : in     Col_Vector;
      Q            : in     Matrix;
      Starting_Col : in     Index;
      Final_Col    : in     Index;
      Max_Err      :    out Real;
      Max_id       :    out Index)
   is
      Err, Max_Eig, Vec_Norm, tst : Real;
      A_times_Eig_r : Col_Vector_e := (others =>  0.0);
      Eigvec_r : Col_Vector_e := (others =>  0.0);
      Err_Vec_r : Col_Vector := (others =>  0.0);
      W_times_EigVec_r : Col_Vector_e := (others =>  0.0);
   begin
      Max_id  := Submatrix_Range'First;
      Max_Err := 0.0;

      Max_Eig := 0.0;
      for Col_id in Submatrix_Range loop
         tst := Abs (W_r(Col_id));
         if tst > Max_Eig then  Max_Eig := tst;  end if;
      end loop;

      for Col_id in Submatrix_Range loop

         for j in Submatrix_Range loop
            Eigvec_r(j) := Real_e (Q(j, Col_id));
         end loop;

         A_times_Eig_r := Product (A, Eigvec_r, Starting_Col, Final_Col);

         for j in Submatrix_Range loop
            W_times_EigVec_r(j) := Real_e(W_r(Col_id)) * Eigvec_r(j);
            Err_Vec_r(j) := Real (A_times_Eig_r(j) - W_times_EigVec_r(j));
         end loop;
      
         -- The vectors Eigvec are already normalized
 
         Vec_Norm := Norm (Err_Vec_r, Starting_Col, Final_Col);

         if Max_Eig > 2.0**(-Real'Machine_Emax) * Vec_Norm  then
            Err := Vec_Norm / (Max_Eig + Min_Allowed_Real);
         else
            Err := 0.0;
         end if;

         -- ||A*v - lambda*v|| / ||lambda_max||  over all normalized v

         if Err > Max_Err then  Max_Err := Err; Max_id := Col_id; end if;

       end loop;

   end Error_in_Decomposition;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7 : string := "") is
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
     new_line;
     begin
        put ("Type a character to continue: ");
        get_immediate (Continue);
     exception
        when others => null;
     end;
   end pause;

   Q      : Matrix;
   A_true : Matrix;
   x, Max_Error : Real := 0.0;
   Max_Id : Index;
begin

   Pause( 
     "Test1: Eigen-decomposition of matrix A.",
     "Error is estimated by measuring Max value of",
     " ",
     "         |A*v - lambda*v| / |max_eigenvalue_of_A|",
     " ",
     "over all normalized eigenvectors v. If 15 digit Reals are used, then",
     "if all goes well, expect the error to be a few parts per 10**15."
     );

   new_line(2);
   put ("Testing QR_Symmetric_Eigen on the following"); 
   put (Integer'Image (Matrix_Size)); put (" x");
   put (Integer'Image (Matrix_Size)); 
   put (" matrices:"); 
   new_line(2);

   for Chosen_Matrix in Matrix_Id loop
 --for Chosen_Matrix in Ding_Dong .. Ding_Dong loop
 --for Chosen_Matrix in All_Ones .. All_Ones loop
 --for Chosen_Matrix in Lower_Ones .. Lower_Ones loop
 --for Chosen_Matrix in vandermonde .. vandermonde loop
 --for Chosen_Matrix in Companion_2 .. Companion_2 loop
 --for Chosen_Matrix in hilbert .. hilbert loop
 --for Chosen_Matrix in random_32_bit .. random_32_bit loop
 --for Chosen_Matrix in Pas_Fib .. Pascal loop
  
      Init_Matrix (A, Chosen_Matrix, Submatrix_Range'First, Submatrix_Range'Last);
  
      --  Want to reuse (unchanged) the error calculation used for Peters_Eigen.
      --  So now we have 2 types of matrices: type Matrix, and type Real_Matrix.
 
--      for c in Submatrix_Range'First .. Submatrix_Range'Last-1 loop
--      for r in c+1 .. Submatrix_Range'Last  loop
--       --A(r, c) := A(c, r);
--         A(c, r) := A(r, c);
--      end loop;
--      end loop;
 

      -- MUST be symmetric

      for c in Submatrix_Range loop
      for r in Submatrix_Range loop
         A_true(r, c) := 0.5 * (A(c, r) + A(r, c));
      end loop;
      end loop;

      A := A_true;

--for i in 1..100 loop
      QR_Method.Eigen_Decompose
        (A                      => A,
         Q                      => Q,
         Eigenvals              => Eigvals,
         Start_Col              => Starting_Col, 
         Final_Col              => Final_Col,
         Eigenvectors_Desired   => True);
--A := A_true;
--end loop;

      x := x + Q(2,2) + A(2,2);

      if true then
         Error_in_Decomposition 
           (A_true, 
            Eigvals,
            Q,
            Starting_Col, Final_Col, 
            Max_Error, Max_id);
  
         put ("Error in decomposition ~"); 
         put (Real'Image(Max_Error));
         put (" for matrix "); put (Matrix_id'Image(Chosen_Matrix));
         new_line;
      end if;
 
      if false then
         put (0.0); put (A(Starting_Col, Starting_Col));
         new_line(1);
         for r in Submatrix_Range'First .. Submatrix_Range'Last-1 loop
            put (A(r+1, r)); put (A(r+1, r+1));
            new_line(1);
         end loop;
      end if;
  
      if false then
         New_Line;
         for I in Submatrix_Range loop
            put (Eigvals(I));
            New_Line(1);
         end loop;
      end if;
  
      if true then
         declare
            Sum : Real_e := +0.0;
          --Product : Real := +1.0;
         begin
            for I in Submatrix_Range loop
               Sum := Sum + Real_e (Eigvals(I));
             --Product := Product * Eigvals(I);
            end loop;
            put ("Eigval sum:"); put (Real (Sum)); 
            New_Line; 
          --put (Product);
         end;
      end if;
  
   end loop;
   put (x);
  
end qr_symmetric_eigen_tst_1;
