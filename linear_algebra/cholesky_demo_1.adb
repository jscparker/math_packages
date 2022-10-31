
-- Test LU decomposition on a
-- real valued square matrix.

with Cholesky_LU;
With Text_IO; use Text_IO;

procedure cholesky_demo_1 is

   type Real is digits 15;

   type Index is range 0..2**7-1;

   Starting_Index : constant Index := Index'First + 0;
   Max_Index      :          Index := Index'Last  - 0;

   type Matrix is array(Index, Index) of Real;
   --  Row major form is appropriate for Matrix * Row_Vector
   --  operations, which dominate the algorithm in procedure Solve.

   package lu is new Cholesky_LU (Real, Index, Matrix);
   use lu;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

   Zero_Vector : constant Row_Vector := (others => 0.0);
   C, C_Inverse : Matrix := (others => (others => 0.0));

   IO_Max_Index : Integer := 4;
   Sum, Max_Error, Del, Error_Off_Diag, Error_Diag : Real;

   type Matrix_Id is (Easy_Matrix, Small_Diagonal, Upper_Ones,
                      Lower_Ones, Moler, Hilbert);

 --Desired_Matrix : Matrix_Id  := Easy_Matrix;
   Desired_Matrix : Matrix_Id  := Small_Diagonal;
 --Desired_Matrix : Matrix_Id  := Upper_Ones;
 --Desired_Matrix : Matrix_Id  := Lower_Ones;
 --Desired_Matrix : Matrix_Id  := Moler;
 --Desired_Matrix : Matrix_Id  := Hilbert;

   type Real_Extended  is digits 15;
 --type Real_Extended  is digits 18;   -- 18 on intel

   e_Sum : Real_Extended;

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8 : string := "") is
     Continue : Character := ' ';
   begin
     new_line;
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

   ------------
   -- Invert --
   ------------
  
   --  Get Inverse of the Matrix:
   
   procedure Invert 
     (M              : in     Matrix;
      M_Inv          :    out Matrix;
      Max_Index      : in     Index;
      Starting_Index : in     Index;
      Max_Error      :    out Real)
   is 
      Solution_Vector : Row_Vector;
      Unit_Vector  : Row_Vector := (others => 0.0);
      Error        : Col_Vector := (others => 0.0);
      M_LU         : Matrix := M;
      Diag_Inverse : Row_Vector;
   begin
      Max_Error := 0.0;
     
      LU_decompose 
        (M_LU,
         Diag_Inverse,
         Max_Index,
         Starting_Index);

      for I in Starting_Index..Max_Index loop

         if I > Starting_Index then 
            Unit_Vector(I-1) := 0.0;
         end if;

         Unit_Vector(I) := 1.0;
         Solve
           (Solution_Vector, 
            Unit_Vector,
            M_LU, 
	    Diag_Inverse,
            Max_Index,
	    Starting_Index);

        --  Solve M*Solution_Vector = Unit_Vector   (for Solution_Vector).

         Error := Unit_Vector - Product (M, Solution_Vector, Max_Index, Starting_Index);
         for I in Starting_Index..Max_Index loop
            if Abs(Error(I)) > Max_Error then
               Max_Error := Abs(Error(I));
            end if;
         end loop;

         -- Solution vector is the I-th column of M_Inverse:

         for J in Starting_Index..Max_Index loop
            M_Inv (J, I) := Solution_Vector(J);
         end loop;

      end loop;
   
   end Invert;
   
begin

   put("Maximum matrix size is "& 
     Integer'Image (Zero_Vector'length-(Integer(Starting_Index)-Integer(Index'First))));
     new_Line;
   put("Input Size Of Matrix To Invert (enter an Integer)"); new_Line;
   get(IO_Max_Index);
   Max_Index := Starting_Index + Index (IO_Max_Index-1);
   
   C := (others => (others => 0.0));

   case Desired_Matrix is
   when Easy_Matrix =>

      for I in Index loop
	 C(I, I) := 1.010101010101;
      end loop;

      for BottomDiagonal in Starting_Index+1..Index'Last loop
      for Row in BottomDiagonal..Index'Last loop
	 C(Row, Row-BottomDiagonal+Starting_Index) 
	         := 0.013 * Real(Row) / Real(Index'Last) 
                           + 0.10101010101 / Real(BottomDiagonal);
      end loop;
      end loop;

      for Row in Starting_Index+1..Index'Last loop
      for Col in Starting_Index..Row-1 loop
         C(Col, Row) :=  C(Row, Col) + 0.333;  
      end loop;
      end loop;

   when Small_Diagonal  =>

      for I in Index loop
	 C(I, I) := 2.0;
      end loop;

      for BottomDiagonal in Starting_Index+1..Index'Last loop
      for Row in BottomDiagonal..Index'Last loop
	 C(Row, Row-BottomDiagonal+Starting_Index) 
	         := 0.013 * Real(Row) / Real(Index'Last) 
                           + 1.0 / Real(BottomDiagonal);
      end loop;
      end loop;

      for Row in Starting_Index+1..Index'Last loop
      for Col in Starting_Index..Row-1 loop
         C(Col, Row) :=  C(Row, Col) + 0.333;  -- this is tough on it.
      end loop;
      end loop;
   
   when Upper_Ones  =>

      C := (others => (others => 0.0));

      for Row in Starting_Index .. Index'Last loop
      for Col in Row .. Index'Last loop
         C(Row, Col) := 1.0; 
      end loop;
      end loop;
   
   when Lower_Ones  =>

      C := (others => (others => 0.0));

      for Row in Starting_Index .. Index'Last loop
      for Col in Row .. Index'Last loop
         C(Col, Row) := 1.0; 
      end loop;
      end loop;
   
   when Moler  =>

      C := (others => (others => 0.0));

      for Row in Starting_Index .. Index'Last loop
      for Col in Starting_Index .. Index'Last loop
         C(Row, Col) := Real(Index'Min(Col,Row)) -  Real(Starting_Index) - 1.0; 
      end loop;
      end loop;
   
      for Col in Starting_Index .. Index'Last loop
	  C(Col, Col) := Real(Col) -  Real(Starting_Index) + 1.0; 
      end loop;
       
   when Hilbert  =>
   
   --  Construct Hilbert's matrix Plus epsilon:
   
      C := (others => (others => 0.0));

      for Row in Starting_Index .. Index'Last loop
      for Col in Starting_Index .. Index'Last loop
	  C(Row, Col) := 
             1.0 / (Real(Row) + Real(Col) - 2.0*Real(Starting_Index) + 1.0);
      end loop;
      end loop;
       
   end case;
   
   --  symmetrize C for Cholesky:

   declare Val : Real; 
   begin
   for Row in Index loop
      if Row < Index'Last then
      for Col in Row+1 .. Index'Last loop
         Val := 0.5 * (C(Col, Row) + C(Row, Col));
         C(Col, Row) := Val;
         C(Row, Col) := Val;
      end loop;
      end if;
   end loop;
   end;

   --  Construct inverse of C:

   Invert  (C, C_Inverse, Max_Index, Starting_Index, Max_Error);

   new_line;
   put("We just took the matrices inverse.  Error estimate follows.");
   new_line;
   put ("Max Error according to Residual function is: "); put (Max_Error);
   new_line;
   
   --  Multiply Original C and C_Inverse as test.  Get Max error:

   Error_Off_Diag := 0.0;
   Error_Diag     := 0.0;
   for I in Starting_Index..Max_Index loop
   for J in Starting_Index..Max_Index loop

      e_Sum := 0.0;
      for K in Starting_Index..Max_Index loop
         e_Sum := e_Sum + Real_Extended (C_Inverse(I, k)) * Real_Extended (C(k, J));
      end loop;
      Sum:= Real(e_Sum);

      --  Product(I,J) := Sum;
      --  The product should be the unit matrix.  Calculate the error:

      if I = J then
         Del := Abs (Sum - 1.0);
         if Del > Error_Diag then
            Error_Diag := Del;
         end if;
      else
         Del := Abs (Sum);
         if Del > Error_Off_Diag then
            Error_Off_Diag := Del;
         end if;
      end if;

   end loop;
   end loop;
   
   pause
     ("We just took the product the original matrix with its inverse,",
      "and then calculated the difference between this product and",
      "the identity matrix.  The difference is printed below. The difference",
      "should be near 10**(-15) only if the matrix is well conditioned.");

   new_line; 
   put("Max difference along diagonals of product:     "); put(Error_Diag);

   new_line; 
   put("Max difference along off-diagonals of product: "); put(Error_Off_Diag);
   new_line;
   
   
   --  Multiply Original C and C_Inverse as test.  Get Max error:

   Error_Off_Diag := 0.0;
   Error_Diag     := 0.0;
   for I in Starting_Index..Max_Index loop
   for J in Starting_Index..Max_Index loop

      e_Sum := 0.0;
      for K in Starting_Index..Max_Index loop
         e_Sum := e_Sum + Real_Extended (C_Inverse(I, k)) * Real_Extended (C(k, J));
      end loop;
      Sum:= Real(e_Sum);

      --  Product(I,J) := Sum;
      --  The product should be the I matrix.  Calculate the error:

      if I = J then
         Del := Abs (Sum - 1.0);
         if Del > Error_Diag then
            Error_Diag := Del;
         end if;
      else
         Del := Abs (Sum);
         if Del > Error_Off_Diag then
            Error_Off_Diag := Del;
         end if;
      end if;

   end loop;
   end loop;
   
   pause
     ("Just took the product of the inverse matrix with the original matrix,",
      "and then calculated the difference between this product and",
      "the identity matrix.  The difference is printed below.",
      "The product of the inverse matrix with the original matrix does not",
      "equal the identity matrix unless the original matrix is well conditioned.");

   new_line; 
   put("The difference along diagonals:     "); put(Error_Diag);

   new_line; 
   put("The difference along off-diagonals: "); put(Error_Off_Diag);
   new_line;
   
end;
