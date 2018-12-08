
---------------------------------------------------------------------------------
-- package body Banded_LU, LU decomposition, equation solving for banded matrices
-- Copyright (C) 1995-2018 Jonathan S. Parker 
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

with Text_IO;

package body Banded_LU is

  -------------
  -- Product --
  -------------
  
  --  Matrix Vector multiplication

  function Product
    (A              : in Banded_Matrix;
     X              : in Column;
     Final_Index    : in Index  := Index'Last;
     Starting_Index : in Index  := Index'First) 
     return Column
  is
     Result : Column := (others => Zero);
     Sum : Real;
     Col_First, Col_Last : Integer;
  begin   
     for Row in Starting_Index .. Final_Index loop
        Sum       := Zero;
        Col_First := Integer'Max (Row - No_Of_Off_Diagonals, Starting_Index);
        Col_Last  := Integer'Min (Row + No_Of_Off_Diagonals, Final_Index);
        for Col in Col_First .. Col_Last loop
           Sum := Sum + A (Row)(Col - Row) * X (Col);
        end loop;
        Result(Row) := Sum;
     end loop;

     return Result;
        
  end Product;
  
  ---------------------
  -- Refine_Solution --
  ---------------------
  
  --  if  No_Of_Iterations=0  then usual solution is returned.   
  --  if  No_Of_Iterations=1  then solution is refined iteratively once.
  --
  --  Not necessarily much use if error is due to ill-conditioning of Matrix A.
  --
  --  Iterative refinement of the solution returned by LU_decompose() and
  --  Solve().  Uses the Newton-like iteration for the solution of A*X = b,
  --
  --     X_{k+1} = X_k + A_Inverse_Approximate * (b - A*X_k).
  --
  --  Here A_Inverse_Approximate (we will call it V below) represents the 
  --  solution returned by LU_decompose() followed by Solve().
  --
  --  if y = exact error in 1st iteration: y = X_inf - X_1, then y is the
  --  exact solution of A*y = d_1 where d_1 = b - A*X_1.
  --  Let V denote approximate inverse of A.  Iterate for y using
  --
  --      Delta_Y_{k+1} == Y_{k+1} - Y_k = V*(d_1 - A*Y_k).
  --
  --  Remember Y = exact error in 1st iteration = SUM (Delta_Y_k's).  
  --  Here's the actual method:
  --
  --   Let d_1 = b - A*X_1 (the standard Residual: 1st estimate of error in A*X = b)
  --   Delta_Y_1 = V*d_1
  --   Let d_2 = d_1 - A*Delta_Y_1
  --   Delta_Y_2 = V*(d_1 - A*Delta_Y_1) = V*d_2 
  --   Let d_3 = d_2 - A*Delta_Y_2
  --   Delta_Y_3 = V*(d_1 - A*Delta_Y_1 - A*Delta_Y_2) = V*d_3
  -- 
  --  so:    d_k = d_{k-1} - A*Delta_Y_{k-1}; Delta_Y_k = V*d_k
  --
  --  Sum the Delta_Y_k's to get the correction to X_1: Y = SUM (Delta_Y_k's).
  --
  procedure Refine_Solution
    (X                :    out Column;
     B                : in     Column;
     A_LU             : in     Banded_Matrix;
     Diag_Inverse     : in     Column;
     A                : in     Banded_Matrix;
     No_Of_Iterations : in     Natural := 1;
     Final_Index      : in     Index   := Index'Last;
     Starting_Index   : in     Index   := Index'First)
  is   
     Delta_Y, X_1, A_Y : Column := (others => 0.0);
     D_k : Column := (others => 0.0);
     Y_f : Column := (others => 0.0);
  begin
     --  Get X_1 as defined above. 

     Solve (X_1, B, A_LU, Diag_Inverse, Final_Index, Starting_Index);
     
     if No_Of_Iterations > 0 then
     
        A_Y := Product(A, X_1, Final_Index, Starting_Index);

        -- D_1:
        for I in Starting_Index .. Final_Index loop
           D_k(I) := B(I) - A_Y(I);
        end loop;

        Solve (Delta_Y, D_k, A_LU, Diag_Inverse, Final_Index, Starting_Index);

        --  Y_f is Sum of all the iterated  Delta_Y's. Initialize it:
        Y_f := Delta_Y;
        
	for Iteration in 1..No_Of_Iterations-1 loop

           --  get d_k = d_k - A*Delta_Y_k
           
           A_Y := Product (A, Delta_Y, Final_Index, Starting_Index);
           for I in Starting_Index .. Final_Index loop
              D_k(I) := D_k(I) - A_Y(I);
           end loop;
 
           --  get Delta_Y = V*D_k:
           
           Solve (Delta_Y, D_k, A_LU, Diag_Inverse, Final_Index, Starting_Index);

           --  Accumulate Y_f: the full correction to X_1:

           for I in Starting_Index .. Final_Index loop
              Y_f(I) := Y_f(I) + Delta_Y(I);
           end loop;

	end loop;

     end if;

     for I in Starting_Index..Final_Index loop
        X(I) := Y_f(I) + X_1(I);
     end loop;
     
  end Refine_Solution;
  
  ----------------
  -- Matrix_Val --
  ----------------

  --  Translates (Row, Col) to (I, Diagonal_id) using
  --  the formula I = Row, and Diagonal_id = Col - Row.
  --
  --  Banded Matrices are by definition 0 everywhere except on the
  --  diagonal bands.  So 0 is returned if (Row, Col) is not in the
  --  banded region.

  function Matrix_Val 
    (A     : Banded_Matrix;
     Row    : Index;
     Col    : Index) 
     return Real
  is
     Diag_ID : constant Integer := (Col - Row);
     Result  : Real;
  begin
     if Abs Diag_ID  > No_Of_Off_Diagonals then
        Result := 0.0;
     else
        Result := A(Row)(Diag_ID);
     end if;
     return Result;
  end;

  ------------------
  -- LU_Decompose --
  ------------------
  
  -- Translates from (Row, Col) indices to (I, Diagonal)
  -- with the formula I = Row, and Diagonal = Col - Row.

  procedure LU_Decompose 
    (A              : in out Banded_Matrix;
     Diag_Inverse   :    out Column;
     Final_Index    : in     Index := Index'Last;
     Starting_Index : in     Index := Index'First)
  is
     Stage : Index;
     Sum : Real;
     Col_First, Row_Last : Integer;
     Min_Allowed_Pivot_Ratio, Min_Allowed_Pivot_Val : Real;
     Reciprocal_Pivot_Val, Pivot_Val, Abs_Pivot_Val : Real;
     Max_Pivot_Val   : Real := Min_Allowed_Real;
     Min_Pivot_Ratio : constant Real := 2.0**(-Real'Machine_Mantissa) * 1.0E-3;
  begin
     Diag_Inverse := (Others => 0.0);
  
     if Final_Index - Starting_Index + 1 < No_Of_Off_Diagonals + 1 then
        text_io.put ("Matrix Size must be >= No_Of_Off_Diagonals+1.");
        raise Constraint_Error;
     end if;
  
     Min_Allowed_Pivot_Ratio := Min_Pivot_Ratio;

    -- Step 0. 1 X 1 matrices: They can't exist because of above.
 
    -- Step 1.  The outer loop. 
    -- At each stage we calculate row "stage" of the Upper matrix U
    -- and Column "Stage" of the Lower matrix L.
    -- The matrix A is overwritten with these, because the elements
    -- of A in those places are never needed in future stages.
    -- However, the elements of L ARE needed in those places,
    -- so to get those elements we will be accessing A (which stores them).
 
    for Stage in Starting_Index..Final_Index-1 loop
    
       Row_Last := Integer'Min (Stage + No_Of_Off_Diagonals, Final_Index);
       
       if Stage > Starting_Index then
       for J in Stage .. Row_Last loop
 
          Sum := 0.0;
          --for K in Starting_Index .. Stage-1 loop
          -- --Sum := Sum + L(J)(K)*U(K)(Stage);
          --   Sum := Sum + A(J)(K)*A(K)(Stage);
          --end loop;
          Col_First := Integer'Max (J - No_Of_Off_Diagonals, Starting_Index);
          for K in Col_First..Stage-1 loop
             Sum := Sum + A(J)(K-J)*A(K)(Stage-K);
          end loop;
 
        --L(J)(Stage)   := L(J)(Stage) - Sum;
        --L(J)(Stage-J) := L(J)(Stage-J) - Sum;
          A(J)(Stage-J) := A(J)(Stage-J) - Sum;
 
       end loop;
       end if;
       
 
       -- Step 2: Get row "stage" of U and
       -- column "stage" of L. Notice these formulas update
       -- only (Stage+1..Last) elements of the respective row
       -- and column, and depend on only (1..Stage) elements
       -- of U and L, which were calculated previously, and stored in A.
 
       Pivot_Val     := A(Stage)(0);
       Abs_Pivot_Val := Abs (Pivot_Val);
 
       if Abs_Pivot_Val > Max_Pivot_Val then
          Max_Pivot_Val := Abs_Pivot_Val;
       end if;
 
       Min_Allowed_Pivot_Val := Max_Pivot_Val*Min_Allowed_Pivot_Ratio + Min_Allowed_Real;
 
       if (Abs_Pivot_Val < Min_Allowed_Pivot_Val) then
          Min_Allowed_Pivot_Val := Real'Copy_Sign (Min_Allowed_Pivot_Val, Pivot_Val);
          Reciprocal_Pivot_Val := 1.0 / Min_Allowed_Pivot_Val;
       else
          Reciprocal_Pivot_Val := 1.0 / Pivot_Val;
       end if;
 
       if (Abs_Pivot_Val < Min_Allowed_Real) then
          Reciprocal_Pivot_Val := 0.0;
       end if;
 
       Diag_Inverse(Stage) := Reciprocal_Pivot_Val;
 
       for J in Stage+1..Row_Last loop
       
          Sum := 0.0;
          if Stage > Starting_Index then
          --for K in Starting_Index .. Stage-1 loop
          --  --Sum := Sum + L(Stage)(K)*U(K)(J);
          --    Sum := Sum + A(Stage)(K)*A(K)(J);
          --end loop;
          Col_First := Integer'Max (Starting_Index, J - No_Of_Off_Diagonals);
          for K in Col_First..Stage-1 loop
             Sum := Sum + A(Stage)(K-Stage) * A(K)(J-K);
          end loop;
          end if;
          
        --U(Stage)(J)       := (A(Stage)(J)       - Sum) * Scale_Factor;
          A(Stage)(J-Stage) := (A(Stage)(J-Stage) - Sum) * Reciprocal_Pivot_Val;

       end loop;
 
     end loop;
 
 
     -- Step 3: Get final row and column.
 
     Stage := Final_Index;
     Sum   := 0.0;
     --for K in Starting_Index .. Stage-1 loop
     --  --Sum := Sum + L(Stage)(K)*U(K)(Stage);
     --    Sum := Sum + A(Stage)(K)*A(K)(Stage);
     --end loop;
     Col_First := Integer'Max(Starting_Index, Integer(Stage)-No_Of_Off_Diagonals);
     for K in Col_First..Stage-1 loop
        Sum := Sum + A(Stage)(K-Stage)*A(K)(Stage-K);
     end loop;
 
     A(Stage)(0) := A(Stage)(0) - Sum;
 
     Pivot_Val     := A(Stage)(0);
     Abs_Pivot_Val := Abs (Pivot_Val);
 
     if Abs_Pivot_Val > Max_Pivot_Val then
        Max_Pivot_Val := Abs_Pivot_Val;
     end if;
 
     Min_Allowed_Pivot_Val := Max_Pivot_Val*Min_Allowed_Pivot_Ratio + Min_Allowed_Real;

     if Abs_Pivot_Val < Min_Allowed_Pivot_Val then
        Min_Allowed_Pivot_Val := Real'Copy_Sign (Min_Allowed_Pivot_Val, Pivot_Val);
        Reciprocal_Pivot_Val := 1.0 / Min_Allowed_Pivot_Val;
     else
        Reciprocal_Pivot_Val := 1.0 / Pivot_Val;
     end if;
 
     if Abs_Pivot_Val < Min_Allowed_Real then
         Reciprocal_Pivot_Val := 0.0;
     end if;

     Diag_Inverse(Stage) := Reciprocal_Pivot_Val;

  end LU_Decompose;

  -----------
  -- Solve --
  -----------

  procedure Solve
    (X              :    out Column;
     B              : in     Column;
     A_LU           : in     Banded_Matrix;
     Diag_Inverse   : in     Column;
     Final_Index    : in     Index := Index'Last;
     Starting_Index : in     Index := Index'First)
  is
     Z             : Column;
     ID_of_1st_non_0 : Index := Starting_Index;
     Sum           : Real;
     Col_First     : Index;
     Col_Last      : Index;
  begin
     for Row in Index loop
        X(Row) := 0.0;
     end loop;
     
     -- An optimization to make matrix inversion efficient. 
     -- in Banded_Matrix inversion, the input vector B is
     -- is a unit vector: it is all zeros except for a 1.0. Need to
     -- to find 1st non-zero element of B:

     for I in Starting_Index..Final_Index loop
        if Abs (B(I)) > 0.0 then
            ID_of_1st_non_0 := I;
            exit;
        end if;
     end loop;

     -- In solving for Z in the equation L Z = B, the Z's will
     -- all be zero up to the 1st non-zero element of B.

     if ID_of_1st_non_0 > Index'First then
        for I in Starting_Index..ID_of_1st_non_0-1 loop
           Z(I) := 0.0;
        end loop;
     end if;

     -- The matrix equation is in the form L * U * X = B.
     -- First assume U * X is Z, and
     -- solve for Z in the equation L Z = B.

     Z(ID_of_1st_non_0) := B(ID_of_1st_non_0) * Diag_Inverse(ID_of_1st_non_0);

     if ID_of_1st_non_0 < Final_Index then
       for Row in ID_of_1st_non_0+1..Final_Index loop
          Sum := 0.0;
          Col_First := Integer'Max (Starting_Index, Row - No_Of_Off_Diagonals);
          for Col in Col_First .. Row-1 loop
              Sum := Sum + A_LU(Row)(Col-Row) * Z(Col);
          end loop;
          Z(Row) := (B(Row) - Sum) * Diag_Inverse(Row);
       end loop;
     end if;

     -- Solve for X in the equation U X = Z.

     X(Final_Index) := Z(Final_Index);
     if Final_Index > Starting_Index then
        for Row in reverse Starting_Index..Final_Index-1 loop
           Sum := 0.0;
           Col_Last := Integer'Min (Final_Index, Row + No_Of_Off_Diagonals);
           for Col in Row+1 .. Col_Last loop
               Sum := Sum + A_LU(Row)(Col-Row) * X(Col);
           end loop;
           X(Row) := (Z(Row) - Sum);
        end loop;
     end if;

  end Solve;

end Banded_LU;
