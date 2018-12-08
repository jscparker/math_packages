
-----------------------------------------------------------------------
-- package body Crout_LU, LU decomposition, with equation solving
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
---------------------------------------------------------------------------------

package body Crout_LU is

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  Min_Allowed_Real : constant Real := Two ** (Real'Machine_Emin - Real'Machine_Emin / 8);

  ---------
  -- "-" --
  ---------
  
  function "-"
    (A, B : in Col_Vector)
     return Col_Vector
  is
     Result : Col_Vector;
  begin   
     for J in Index loop
        Result(J) := A(J) - B(J);  
     end loop;
     return Result;
  end "-";
  
  -------------
  -- Product --
  -------------

  function Product
    (A              : Matrix;
     V              : Row_Vector;
     Final_Index    : Index := Index'Last;
     Starting_Index : Index := Index'First)
     return Row_Vector
  is
     Result : Col_Vector := (others => Zero);
     Sum : Real;
  begin
    for i in Starting_Index .. Final_Index loop
       Sum := Zero;
       for j in Starting_Index .. Final_Index loop
          Sum := Sum + A(i, j) * V(j);
       end loop;
       Result (i) := Sum;
    end loop;

    return Result;

  end Product;

  --------------------------
  -- Scale_Cols_Then_Rows --
  --------------------------

  procedure Scale_Cols_Then_Rows
    (A              : in out Matrix;
     Scalings       :    out Scale_Vectors;
     Final_Index    : in     Index := Index'Last;
     Starting_Index : in     Index := Index'First) 
  is
     Sum, Scale_Factor : Real;
     Power_of_Two : Integer;
  begin
     --  Scale each column to near unity:

     Scalings (For_Cols) := (others => One);

     for Col in Starting_Index .. Final_Index loop
        Sum := Zero; 
        for j in Starting_Index .. Final_Index loop
           Sum := Sum + Abs A(j, Col);  
        end loop;

        Power_of_Two := Real'Exponent (Sum + Min_Allowed_Real);
        Scale_Factor := Two ** (-Power_of_Two);

        for j in Starting_Index .. Final_Index loop
           A(j, Col) := Scale_Factor * A(j, Col);  
        end loop;

        Scalings (For_Cols)(Col) := Scale_Factor;
     end loop;

     --  Scale each row to near unity:

     Scalings (For_Rows) := (others => One);

     for Row in Starting_Index .. Final_Index loop
        Sum := Zero; 
        for j in Starting_Index .. Final_Index loop
           Sum := Sum + Abs A(Row, j);  
        end loop;

        Power_of_Two := Real'Exponent (Sum + Min_Allowed_Real);
        Scale_Factor := Two ** (-Power_of_Two);

        for j in Starting_Index .. Final_Index loop
           A(Row, j) := Scale_Factor * A(Row, j);  
        end loop;

        Scalings (For_Rows)(Row) := Scale_Factor;
     end loop;

  end Scale_Cols_Then_Rows;
  
  ------------------
  -- LU_Decompose --
  ------------------

  -- The upper matrix is U, the lower L.
  -- We assume that the diagonal elements of L are One. Thus
  -- the diagonal elements of U but not L appear on the
  -- the diagonal of the output matrix A.

  procedure LU_Decompose
    (A                : in out Matrix;
     Scalings         :    out Scale_Vectors;
     Row_Permutation  :    out Rearrangement;
     Final_Index      : in     Index   := Index'Last;
     Starting_Index   : in     Index   := Index'First;
     Scaling_Desired  : in     Boolean := False)
  is
     Stage : Index;
     tmp_Index, The_Pivotal_Row : Index;
     Sum, tmp : Real;

     Min_Allowed_Pivot_Val, Reciprocal_Pivot_Val : Real;
     Pivot_Val, Abs_Pivot_Val : Real;

     Min_Pivot_Ratio : constant Real := Two**(-Real'Machine_Mantissa-24);

     Max_Pivot_Val : Real := Min_Allowed_Real;

     -----------------------------
     -- Find_Max_Element_Of_Col --
     -----------------------------

     procedure Find_Max_Element_Of_Col
       (Col_ID                  : in Index;
        Starting_Index          : in Index;
        Index_of_Max_Element    : out Index;
        Val_of_Max_Element      : out Real;
        Abs_Val_of_Max_Element  : out Real)
     is
        Pivot_Val, Abs_Pivot_Val : Real;
     begin
        Val_of_Max_Element     := A (Starting_Index, Col_ID);
        Abs_Val_of_Max_Element := Abs (Val_of_Max_Element);
        Index_of_Max_Element   := Starting_Index;
        if Final_Index > Starting_Index then
        for k in Starting_Index+1..Final_Index loop
           Pivot_Val     := A (k, Col_ID);
           Abs_Pivot_Val := Abs (Pivot_Val);
           if Abs_Pivot_Val > Abs_Val_of_Max_Element then
              Val_of_Max_Element     := Pivot_Val;
              Abs_Val_of_Max_Element := Abs_Pivot_Val;
              Index_of_Max_Element   := k;
           end if;
        end loop;
        end if;
     end Find_Max_Element_Of_Col;

  begin

     for I in Index loop
        Row_Permutation(I) := I;
     end loop;

     Scalings(Diag_Inverse) := (others => Zero);
     Scalings(For_Cols)     := (others => One);
     Scalings(For_Rows)     := (others => One);

     if Scaling_Desired then
        Scale_Cols_Then_Rows (A, Scalings, Final_Index, Starting_Index);
     end if;

     -- Step 0: 1 X 1 matrices:

     if Final_Index = Starting_Index then
        Pivot_Val := A(Starting_Index, Starting_Index);
        if Abs (Pivot_Val) <  Min_Allowed_Real then
           A(Starting_Index, Starting_Index) := Zero;
        else
           A(Starting_Index, Starting_Index)      := Pivot_Val;
           Scalings(Diag_Inverse)(Starting_Index) := One / Pivot_Val;
        end if;
        return;
     end if;

     -- Process goes through stages Starting_Index..Final_Index.
     -- The last stage is a special case.
     --
     -- At each stage calculate row "stage" of the Upper
     -- matrix U and Col "Stage" of the Lower matrix L.
     -- The matrix A is overwritten with these, because the elements
     -- of A in those places are never needed in future stages.
     -- However, the elements of U and L ARE needed in those places,
     -- so to get those elements we access A (which stores them).

     for Stage in Starting_Index .. Final_Index-1 loop

        if Stage > Starting_Index then
        for Row in Stage .. Final_Index loop
           Sum := Zero;
           for K in Starting_Index .. Stage-1 loop
            --Sum := Sum + L(Row, K)*U(K, Stage);
              Sum := Sum + A(Row, K)*A(K, Stage);
           end loop;
           A(Row, Stage) := A(Row, Stage) - Sum;
        end loop;
        end if;

        -- Step 2. Swap rows of L and A if necessary.
        -- Do it by swapping rows of A.
        -- Notice that the Rows of U that have already been calculated and
        -- stored in A, namely (1..Stage-1), are untouched by the swap.

        Find_Max_Element_Of_Col
          (Col_ID                 => Stage,
           Starting_Index         => Stage,
           Index_of_Max_Element   => The_Pivotal_Row,
           Val_of_Max_Element     => Pivot_Val,
           Abs_Val_of_Max_Element => Abs_Pivot_Val);

        if The_Pivotal_Row /= Stage then
           for j in Starting_Index .. Final_Index loop
             tmp                   := A(The_Pivotal_Row, j);
             A(The_Pivotal_Row, j) := A(Stage, j);
             A(Stage, j)           := tmp;
	   end loop;

           tmp_Index                        := Row_Permutation(The_Pivotal_Row);
           Row_Permutation(The_Pivotal_Row) := Row_Permutation(Stage);
           Row_Permutation(Stage)           := tmp_Index;
        end if;

        -- Step 3:
        -- Update Ith_row = Stage of the upper triangular matrix U.
        -- Update Ith_col = Stage of the lower triangular matrix L.
        -- The rules are that the diagonal elements of L are 1 even
        -- though Pivot_Val * Reciprocal_Pivot_Val /= 1.
        -- Constraint is that L*U = A when possible.

        if Abs_Pivot_Val > Max_Pivot_Val then
           Max_Pivot_Val := Abs_Pivot_Val;
        end if;

        Min_Allowed_Pivot_Val := Max_Pivot_Val * Min_Pivot_Ratio + Min_Allowed_Real;

        if (Abs_Pivot_Val < Abs Min_Allowed_Pivot_Val) then
           Reciprocal_Pivot_Val  := Zero;
        else
           Reciprocal_Pivot_Val := One / Pivot_Val;
        end if;

        Scalings(Diag_Inverse)(Stage) := Reciprocal_Pivot_Val;
        A(Stage, Stage)               := Pivot_Val;

        for Row in Stage+1 .. Final_Index loop
           A(Row, Stage) := A(Row, Stage) * Reciprocal_Pivot_Val;
        end loop;

        if Stage > Starting_Index then
        for Col in Stage+1 .. Final_Index loop
           Sum := Zero;
           for K in Starting_Index .. Stage-1 loop
             --Sum := Sum + L(Stage, K)*U(K, Col);
               Sum := Sum + A(Stage, K)*A(K, Col);
           end loop;
         --U(Stage, Col) := A(Stage, Col) - Sum;
           A(Stage, Col) := A(Stage, Col) - Sum;
        end loop;
        end if;

     end loop; -- Stage

     -- Step 4: Get final row and column.

     Stage := Final_Index;

     Sum := Zero;
     for K in Starting_Index .. Stage-1 loop
      --Sum := Sum + L(Stage, K)*U(K, Stage);
        Sum := Sum + A(Stage, K)*A(K, Stage);
     end loop;

     Pivot_Val     := A(Stage, Stage) - Sum;
     Abs_Pivot_Val := Abs Pivot_Val;

     Min_Allowed_Pivot_Val := Max_Pivot_Val * Min_Pivot_Ratio + Min_Allowed_Real;

     if (Abs_Pivot_Val < Abs Min_Allowed_Pivot_Val) then
        Reciprocal_Pivot_Val  := Zero;
     else
        Reciprocal_Pivot_Val := One / Pivot_Val;
     end if;

     Scalings(Diag_Inverse)(Stage) := Reciprocal_Pivot_Val;
     A(Stage, Stage)               := Pivot_Val;

  end LU_Decompose;

  --------------
  -- LU_Solve --
  --------------

  procedure LU_Solve
    (X                :    out Row_Vector;
     B                : in     Row_Vector;
     A_LU             : in     Matrix;
     Scalings         : in     Scale_Vectors;
     Row_Permutation  : in     Rearrangement;
     Final_Index      : in     Index := Index'Last;
     Starting_Index   : in     Index := Index'First)
  is
     Z, B2, B_s  : Row_Vector;
     Sum : Real;
  begin
     X := (others => Zero);

     --  A*X = B was changed to (S_r*A*S_c) * (S_c^(-1)*X) = (S_r*B).
     --  The matrix LU'd was (S_r*A*S_c). Let B_s = (S_r*B). Solve for
     --  X_s = (S_c^(-1)*X).
     --
     --  The matrix equation is now  P*L*U*X_s = B_s (the scaled B).
     --
     --  First assume L*U*X_s is B2, and solve for B2 in the equation P*B2 = B_s.
     --  Permute the elements of the vector B_s according to "Permutation":

     -- Get B_s = S_r*B:

     for i in Starting_Index .. Final_Index loop
        B_s(i) := Scalings(For_Rows)(i) * B(i);
     end loop;

     -- Get B2 by solving P*B2 = B_s = B:

     for i in Starting_Index .. Final_Index loop
        B2(i) := B_s(Row_Permutation(i));
     end loop;

     --  The matrix equation is now  L*U*X = B2.
     --  Assume U*X is Z, and solve for Z in the equation L*Z = B2.
     --  Remember that by assumption L(I, I) = One, U(I, I) /= One.

     Z(Starting_Index) := B2(Starting_Index);
 
     if Starting_Index < Final_Index then
	for Row in Starting_Index+1 .. Final_Index loop
           Sum := Zero;
	   for Col in Starting_Index .. Row-1 loop
              Sum := Sum + A_LU(Row, Col) * Z(Col);
	   end loop;
           Z(Row) := B2(Row) - Sum;
	end loop;
     end if;

     --  Solve for X_s in the equation U X_s = Z.

     X(Final_Index) := Z(Final_Index) * Scalings(Diag_Inverse)(Final_Index);

     if Final_Index > Starting_Index then
	for Row in reverse Starting_Index .. Final_Index-1 loop
           Sum := Zero;
	   for Col in Row+1 .. Final_Index loop
              Sum := Sum + A_LU(Row, Col) * X(Col);
	   end loop;
           X(Row) := (Z(Row) - Sum) * Scalings(Diag_Inverse)(Row);
	end loop;
     end if;
     
     --  Solved for the scaled X_s (but called it X); now get the real X:

     for i in Starting_Index .. Final_Index loop
        X(i) := Scalings(For_Cols)(i) * X(i);
     end loop;

  end LU_Solve;

end Crout_LU;
