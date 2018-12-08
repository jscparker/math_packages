
-- Test LU decomposition on a real valued square matrix.

with Ada.Numerics.Generic_elementary_functions;
with Crout_LU;
with Text_IO; use Text_IO;
with Test_Matrices;

procedure crout_lu_tst_2 is

   type Real is digits 15;

   subtype Index is Integer range 0..50;

   Starting_Index : constant Index := Index'First + 0;
   Final_Index    :          Index := Index'Last  - 0;

   type Matrix is array(Index, Index) of Real;

   package math is new Ada.Numerics.Generic_elementary_functions(Real);
   use math; --for sqrt
   package lu is new crout_lu (Real, Index, Matrix);
   use lu;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;
   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;


   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;
   Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin / 4);

   Zero_Vector : constant Row_Vector := (others => Zero);

   A, B, Err : Matrix := (others => (others => Zero));
   Identity, B_inv : Matrix := (others => (others => Zero));
   Relative_Err, Max_Error   : Real;
   IO_Final_Index : Integer := 4;

   Scale_the_Matrix : constant Boolean := True;

   Scale  : Scale_Vectors;

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8 : string := "") is
     Continue : Character := ' ';
   begin
     New_Line;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     if S8 /= "" then put_line (S8); end if;
     new_line;
     begin
        put ("Enter a character to continue: ");
        get_immediate (Continue);
        new_line;
     exception
	when others => null;
     end;
   end pause;

   -----------------------------------
   -- Transpose_of_Left_Times_Right --
   -----------------------------------
  
   function Transpose_of_Left_Times_Right
     (A, B         : in     Matrix;
      Final_Row    : in     Index := Final_Index; 
      Final_Col    : in     Index := Final_Index;
      Starting_Row : in     Index := Starting_Index; 
      Starting_Col : in     Index := Starting_Index)
      return Matrix
   is
      Sum : Real := Zero;
      Result : Matrix := (others => (others => Zero));
   begin
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         Sum := Zero;
         for k in Starting_Col .. Final_Col  loop
            Sum := Sum + A(k, Row) * B(k, Col);
         end loop;
	 Result(Row, Col) := Sum;
      end loop;
      end loop;
      return Result;
   end Transpose_of_Left_Times_Right;

   -------------
   -- Product --
   -------------
  
   function Product
     (A, B         : in     Matrix;
      Final_Row    : in     Index := Final_Index; 
      Final_Col    : in     Index := Final_Index;
      Starting_Row : in     Index := Starting_Index; 
      Starting_Col : in     Index := Starting_Index)
      return Matrix
   is
      Sum : Real := Zero;
      Result : Matrix := (others => (others => Zero));
   begin
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         Sum := Zero;
         for k in Starting_Col .. Final_Col  loop
            Sum := Sum + A(Row, k) * B(k, Col);
         end loop;
	 Result(Row, Col) := Sum;
      end loop;
      end loop;
      return Result;
   end Product;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix;
      Final_Row    : in     Index := Final_Index; 
      Final_Col    : in     Index := Final_Index;
      Starting_Row : in     Index := Starting_Index; 
      Starting_Col : in     Index := Starting_Index)
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

      Scaling := One / (Max_A_Val + Min_Allowed_Real);

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
   
   procedure Invert 
     (M              : in     Matrix;
      M_Inv          :    out Matrix;
      Max_Error      :    out Real)
    --Final_Index    : in     Index;
    --Starting_Index : in     Index)
   is 
      Solution_Vector : Row_Vector;
      Unit_Vector  : Row_Vector := (others => 0.0);
      Error        : Col_Vector := (others => 0.0);
      Scale        : Scale_Vectors;
      Permute      : Rearrangement;
      M_LU         : Matrix := M;
   begin
      Max_Error := 0.0;
     
      LU_decompose 
        (A               => M_LU, 
         Scalings        => Scale, 
         Row_Permutation => Permute,
         Final_Index     => Final_Index,
         Starting_Index  => Starting_Index,
         Scaling_Desired => Scale_the_Matrix);
   
      for I in Starting_Index..Final_Index loop

         if I > Starting_Index then 
            Unit_Vector(I-1) := 0.0;
         end if;

         Unit_Vector(I) := 1.0;

         --  Solve  A*X = B:

         LU_Solve
           (X               => Solution_Vector, 
            B               => Unit_Vector,
	    A_LU            => M_LU, 
            Scalings        => Scale,
	    Row_Permutation => Permute,
            Final_Index     => Final_Index,
	    Starting_Index  => Starting_Index);

        --  Solve M*Solution_Vector = Unit_Vector   (for Solution_Vector).

         Error := Unit_Vector - Product (M, Solution_Vector,Final_Index,Starting_Index);
         for I in Starting_Index..Final_Index loop
            if Abs(Error(I)) > Max_Error then
               Max_Error := Abs(Error(I));
            end if;
         end loop;

         -- Solution vector is the I-th column of M_Inverse:

         for J in Starting_Index..Final_Index loop
            M_Inv (J, I) := Solution_Vector(J);
         end loop;

      end loop;
   
   end Invert;
   
begin

   for Col in Index loop
      Identity(Col, Col) := One;
   end loop;

   put("Maximum matrix size is "& 
      Integer'Image (Zero_Vector'length-(Integer(Starting_Index)-Integer(Index'First))));
      new_Line;
   put("Input Size Of Matrix To Invert (enter an Integer)"); new_Line;
   get(IO_Final_Index);
   Final_Index := Starting_Index + Index (IO_Final_Index-1);
   
   Pause( 
     "Test 1: LU Decomposition of matrix A = P*L*U is usually successful, but if A is",
     "singular of ill-conditioned then the P*L*U is likely to be useless.",
     "If A is singular or ill-conditioned then it's better to solve the normal equations,",
     "in other words solve A'*A * X = A'*b instead of A * X = b, (where A' = Transpose(A)).",
     "A'*A may still be ill-conditioned, but adding 'epsilon*I' to A'*A reliably removes the",
     "singularity. Sometimes iterative refinement is used to get better solutions to A*X = b,",
     "but in this test no iterative refinement is performed."
     );
 
   Pause( 
     "Final remark: if you work with the normal equations, A'*A * X = A'*b, then",
     "Choleski decomposition is more appropriate, but the idea here is to test the",
     "LU decomposition, LU_solve.",
     "The matrices used in the following test are mostly ill-conditioned, so in most cases",
     "the error does not approach machine precision (1.0e-16).",
     " ",
     "Equation solving is done with matrix B = A'*A + epsilon*I."
     );
   new_line;
 
   for Chosen_Matrix in Matrix_id loop
   
      put("For matrix originally of type "); 
      put(Matrix_id'Image(Chosen_Matrix)); Put(":"); 
      new_line; 
 
      -- Get a non-singular B = A'*A + epsilon*I:
 
      Init_Matrix (A, Chosen_Matrix, Starting_Index, Final_Index);
    
      Scale_Cols_Then_Rows
        (A              => A,
         Scalings       => Scale,
         Final_Index    => Final_Index,
         Starting_Index => Starting_Index);
 
      B := Transpose_of_Left_Times_Right (A, A); -- get B = A'*A
 
      --  Now B's eig vals are all >= 0.  Add Eps*Identity to B:  this
      --  shifts all eigs plus Eps = Two**(-28), so all eigs >= Two**(-28):
 
      for Col in Index loop
         B(Col, Col) :=  Two**(-20) + B(Col, Col);
      end loop;
 
      -- Get A_inverse
   
      Invert (B, B_inv, Max_Error);
    
      for Row in Starting_Index .. Final_Index loop
      for Col in Starting_Index .. Final_Index loop
         Err (Row, Col) := Identity(Row, Col) - Product (B, B_inv)(Row, Col);
      end loop;
      end loop;
 
      --Err := Transpose_of_Left_Times_Right (A, Err); 
 
      Relative_Err := Frobenius_Norm (Err) / (Frobenius_Norm (Identity));
   
      put(" Err in I - B*B_inverse is ~ ||I - B*B_inverse|| / ||I|| ="); 
      put(Relative_Err);
      new_line; 
   
      for Row in Starting_Index .. Final_Index loop
      for Col in Starting_Index .. Final_Index loop
         Err (Row, Col) := Identity(Row, Col) - Product (B_inv, B)(Row, Col);
      end loop;
      end loop;
   
      Relative_Err := Frobenius_Norm (Err) / (Frobenius_Norm (Identity));
   
      put(" Err in I - B_inverse*B is ~ ||I - B_inverse*B|| / ||I|| ="); 
      put(Relative_Err);
      new_line; 
   
   end loop; -- Matrix_id
      
end;
