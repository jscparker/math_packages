
-- Tests LU decomposition on a real valued banded matrix.
-- 
-- Shows how interative refinement can clean up the mess when no
-- pivoting is used in LU decomposition.  Also shows that iterative 
-- refinement is ineffective if the problem is ill-conditioned.
-- (Pivoting doesn't solve the problem of ill-conditioning either.)

with Banded_LU;
With Text_IO; use Text_IO;

procedure banded_lu_demo_1 is

   type Real is digits 15;

   -- Create a data structure for making banded matrices.

   No_Of_Off_Diagonals : constant Integer := 20; -- if 1, it's tri-diagonal
   Max_Size_Of_Matrix  : constant Integer := 150;

   package lu is new Banded_LU (Real, Max_Size_Of_Matrix, No_Of_Off_Diagonals);
   use lu;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

   type Matrix is array(Index) of Column; -- not type Banded_Matrix

   Starting_Index : constant Index := Index'First + 10;
   Max_Index      : constant Index := Index'Last  - 10;

   C_band    : Banded_Matrix := (others => (others => 0.0)); -- exported by LU
   C_full    : Matrix        := (others => (others => 0.0));
   C_Inverse : Matrix        := (others => (others => 0.0));

   No_Of_Iterations : Integer := 0;
   Off_Factor      : Real;
   Sum, Del, Error_Off_Diag, Error_Diag, Max_Error : Real;

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8,s9 : string := "") is
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
     if S9 /= "" then put_line (S9); end if;
     new_line(1);
     begin
        put ("Type a character to continue: ");
        get_Immediate (Continue);
        new_line(1);
     exception
        when others => null;
     end;
   end pause;

   ------------
   -- Invert --
   ------------

   --  Get Inverse of the Matrix:

   procedure Invert
     (M                : in     Banded_Matrix;
      M_Inv            :    out Matrix;
      Max_Index        : in     Index;
      Starting_Index   : in     Index;
      Max_Error        :    out Real;
      No_Of_Iterations : in     Integer)
   is
      Unit_Vector, Solution_Vector, Error, Diag_Inverse : Column;
      Zero_Vector : constant Column := (others => 0.0);
      M_LU : Banded_Matrix := M;
   begin
      Max_Error := 0.0;

      LU_decompose (M_LU, Diag_Inverse, Max_Index, Starting_Index);

      for I in Starting_Index..Max_Index loop

          Unit_Vector    := Zero_Vector;
          Unit_Vector(I) := 1.0;
        --Solve(Solution_Vector, Unit_Vector, M_LU, Diag_Inverse, M, 
	--                                    Max_Index, Starting_Index);

          Refine_Solution
             (Solution_Vector, Unit_Vector, M_LU, Diag_Inverse, M,
                             No_Of_Iterations, Max_Index, Starting_Index);

          Error := Product (M, Solution_Vector, Max_Index, Starting_Index);

          for J in Starting_Index .. Max_Index loop
             Error (J) := Error(J) - Unit_Vector(J);
          end loop;

          for I in Starting_Index .. Max_Index loop
             if Abs(Error(I)) > Max_Error then
                Max_Error := Abs(Error(I));
             end if;
          end loop;

          -- Solution vector is the I-th column of M_Inverse:

          for J in Starting_Index .. Max_Index loop
             M_Inv (J)(I) := Solution_Vector(J);
          end loop;

      end loop;

   end Invert;

   -------------------------------------
   -- Test_Inversion_On_Banded_Matrix --
   -------------------------------------

   procedure Test_Inversion_On_Banded_Matrix
     (Off_Factor       : in Real;
      No_Of_Iterations : in Integer)
   is
   begin
      -- construct matrix.

      C_band := (others => (others => 0.0));
      C_full := (others => (others => 0.0));

      for I in Index loop
         C_band(I)(0) := 1.01010101010;
         C_full(I)(I) := C_band(I)(0);
      end loop;

      for Diagonal in Diag_Index'First .. -1 loop
      for Row in Index'First - Diagonal .. Index'Last loop
         C_full(Row)(Diagonal + Row)
                 := 0.033 * Real(Row) / Real(Index'Last)
                               + Off_Factor / Abs Real(Diagonal);
         C_band(Row)(Diagonal)
                 := 0.033 * Real(Row) / Real(Index'Last)
                               + Off_Factor / Abs Real(Diagonal);
      end loop;
      end loop;

      for Diagonal in 1 .. Diag_Index'Last loop
      for Row in Index'First .. Index'Last - Diagonal loop
         C_full(Row)(Diagonal + Row)
                 := 0.013 * Real(Row) / Real(Index'Last)
                               + Off_Factor / Abs Real(Diagonal) + 0.3;
         C_band(Row)(Diagonal)
                 := 0.013 * Real(Row) / Real(Index'Last)
                               + Off_Factor / Abs Real(Diagonal) + 0.3;
      end loop;
      end loop;

      --  Test Matrix_Val procedure:

      for Row in Index loop
      for Col in Index loop
         if Abs (C_full(Row)(Col) - Matrix_Val (C_band, Row, Col)) > 1.0e14 then
            put_line("Some error in procedure Matrix_Val or ... ");
         end if;
      end loop;
      end loop;

     --  Construct inverse of C_band:

     Invert (C_band, C_Inverse, Max_Index, Starting_Index, Max_Error, No_Of_Iterations);

     pause("We just took the matrix inverse.  Error estimate follows.");

     new_line(2);
     put ("Max Error according to Residual function is: "); put (Max_Error);
     new_line;

     --  Multiply Original C_band and C_Inverse as test.  Get Max error:

     Error_Off_Diag := 0.0;
     Error_Diag     := 0.0;
     for I in Starting_Index..Max_Index loop
     for J in Starting_Index..Max_Index loop

        Sum := 0.0;
        for K in Starting_Index..Max_Index loop
           Sum := Sum + C_Inverse(I)(k) * C_full(k)(J);
        end loop;

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
       ("We just took the product the inverse matrix with the original,",
        "matrix and then calculated the difference between this product",
        "and the unit matrix.  The error is printed below.");

     new_line;
     put("Max error along diagonals of product:     "); put(Error_Diag);

     new_line;
     put("Max error along off-diagonals of product: "); put(Error_Off_Diag);
     new_line;

     --  Multiply Original C and C_Inverse as test.  Get Max error:

     Error_Off_Diag := 0.0;
     Error_Diag     := 0.0;
     for I in Starting_Index..Max_Index loop
     for J in Starting_Index..Max_Index loop

        Sum := 0.0;
        for K in Starting_Index..Max_Index loop
           Sum := Sum + C_full(I)(K) * C_Inverse(K)(J);
        end loop;

        --  Product(I,J) := Sum;
        --  The product should be the unit matrix.  Calculate the error:

        if I = J then
           Del := Abs (Sum - 1.0);
           if Del > Error_Diag then
              Error_Diag := Del;
           end if;
        ELSE
           Del := Abs (Sum);
           if Del > Error_Off_Diag then
              Error_Off_Diag := Del;
           end if;
        end if;

     end loop;
     end loop;

     pause
       ("We just took the product the original matrix with the inverse matrix,",
        "and then calculated the difference between this product and",
        "the unit matrix.  The error is printed below.");

     new_line;
     put("Max error along diagonals of product:     "); put(Error_Diag);

     new_line;
     put("Max error along off-diagonals of product: "); put(Error_Off_Diag);
     new_line;

  end Test_Inversion_On_Banded_Matrix;

begin

   --put("Input Size Of Matrix To Invert (enter an Integer)"); new_line;
   --get(IO_Max_Index);
   --Max_Index := Index'First + Index (IO_Max_Index-1);


   pause
     ("The first example demonstrates LU decomposition and inversion of",
      "of a diagonally dominant matrix.   The off-diagonals of the matrix",
      "are all multiplied by 0.5 to make them fall off fast enough for",
      "diagonal dominance (i.e. Off_Factor = 0.5).  No iterative refinement",
      "is done.  (i.e. No_Of_Iterations = 0).");

   Off_Factor       := 0.5;
   No_Of_Iterations := 0;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);

   pause
     ("The 2nd example demonstrates LU decomposition and inversion of",
      "a non-diagonally dominant matrix.   The off-diagonals of the matrix",
      "are all multiplied by 10_000_000_000.0 to make them far greater in value",
      "than the diagonal (i.e. Off_Factor = 1.0e10).  No iterative refinement",
      "is done.  (i.e. No_Of_Iterations = 0).  A large amount of error should",
      "be evident.");

   Off_Factor       := 10_000_000_000.0;
   No_Of_Iterations := 0;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);

   pause
     ("The 3rd example demonstrates LU decomposition and inversion of",
      "a non-diagonally dominant matrix along with iterative refinement. The",
      "off-diagonals of the matrix are again multiplied by 1.0e10 to make",
      "them far greater in value than the diagonal.",
      "Two iterations of the iterative refinement procedure are performed. ",
      "(i.e. No_Of_Iterations = 2).  Error should be greatly reduced.");

   Off_Factor       := 10_000_000_000.0;
   No_Of_Iterations := 2;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);

   pause
     ("The 4rth example demonstrates LU decomposition and inversion of",
      "of a non-diagonally dominant matrix which is mildly ill-conditioned.",
      "The off-diagonals of the matrix are multiplied by 0.8.",
      "No iterative refinement is performed. (i.e. No_Of_Iterations = 0).",
      "Large error should be evident.");

   Off_Factor       := 0.8;
   No_Of_Iterations := 0;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);


   pause
     ("The final example demonstrates LU decomposition and inversion of",
      "of a mildly ill-conditioned matrix along with iterative refinement.",
      "The off-diagonals of the matrix are again multiplied by 0.8.",
      "Iterative refinement is performed. (i.e. No_Of_Iterations = 1).",
      "No amount of iterative refinement helps much here.");

   Off_Factor       := 0.8;
   No_Of_Iterations := 1;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);


   new_line(2);
   put ("Choose your own values for Off_Factor and No_Of_Iterations.");
   new_line(2);
   put ("Input factor for scaling the off diagonals (e.g. 1.0e13):");
   get(Off_Factor);

   new_line;
   put("Input number of iterative refinements of the solution (e.g. 4):");
   get(No_Of_Iterations);
   new_line;

   Test_Inversion_On_Banded_Matrix (Off_Factor, No_Of_Iterations);

end;
