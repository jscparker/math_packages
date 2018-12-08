
--------------------------------------------------------------------------------
-- package body Cholesky_LU, Cholesky decomposition
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

with Ada.Numerics.Generic_Elementary_Functions;

package body Cholesky_LU is

  package Math is new Ada.Numerics.Generic_Elementary_Functions(Real);-- for Sqrt
  use Math;

  --------------
  -- Product --
  --------------
  
  --  Matrix Vector multiplication

  function Product
    (A              : in Matrix;
     X              : in Col_Vector;
     Final_Index    : in Index := Index'Last;
     Starting_Index : in Index := Index'First) 
     return Col_Vector
  is
     Result : Col_Vector := (others => Zero);
     Sum : Real;
  begin   
     for i in Starting_Index .. Final_Index loop 
        Sum := Zero; 
        for j in Starting_Index .. Final_Index loop
           Sum := Sum + A(i, j) * X(j);  
        end loop;
        Result(I) := Sum;
     end loop;

     return Result;
        
  end Product;
  
  ---------
  -- "-" --
  ---------
  
  function "-"(A, B : in Col_Vector) return Col_Vector is
     Result : Col_Vector;
  begin   
     for J in Index loop
        Result(J) := A(J) - B(J);  
     end loop;
     return Result;
  end "-";
  
  ------------------
  -- LU_Decompose --
  ------------------

  procedure LU_Decompose 
    (A               : in out Matrix;
     Diag_Inverse    :    out Col_Vector;
     Final_Index     : in     Index   := Index'Last;
     Starting_Index  : in     Index   := Index'First)
  is
     Final_Stage  : Index;
     Sum, Pivot_Info, Scale_Factor : Real;
  begin

     for I in Starting_Index .. Final_Index loop
        if A(I, I) <= Zero then
           raise Constraint_Error with "Matrix doesn't seem to be positive definite.";
        end if;
     end loop;
     --  Positive definite matrices are never zero or negative on the
     --  diagonal.  The above catches some common input errors, if you want.

     -- 1 X 1 matrix:

     if Final_Index = Starting_Index then
        A(Final_Index, Final_Index) := Sqrt (A(Final_Index, Final_Index));
        Diag_Inverse(Final_Index)   := One / A(Final_Index, Final_Index);
        return;
     end if;

     -- Perform decomp in stages Starting_Index..Final_Index.
     -- The last stage, Final_Index, is a special case.
     -- At each stage calculate row "stage" of the Upper matrix
     -- U and Col "Stage" of the Lower matrix L.

     for Stage in Starting_Index .. Final_Index-1 loop

        Sum := Zero;
        if Stage > Starting_Index then
        for Col in Starting_Index .. Stage-1 loop
          --Sum := Sum + L(Stage, Col) * U(Col, Stage);-- use U = transpose(L)
            Sum := Sum + A(Stage, Col) * A(Stage, Col);
        end loop;
        end if;
        Pivot_Info := A(Stage, Stage) - Sum;

        if Pivot_Info < Min_Allowed_Real then
           raise Constraint_Error with "Matrix doesn't seem to be positive definite.";
        end if;
        
        A(Stage, Stage)     := Sqrt (Pivot_Info);
        Scale_Factor        := One / A(Stage, Stage);
        Diag_Inverse(Stage) := Scale_Factor;

        if Stage > Starting_Index then
           for Row in Stage+1..Final_Index loop
              Sum := Zero;
              for Col in Starting_Index .. Stage-1 loop
               --Sum := Sum + L(Stage, Col) * U(Col, Row); -- use U = transpose(L)
                 Sum := Sum + A(Stage, Col) * A(Row, Col);
              end loop;
              A(Row, Stage) := Scale_Factor * (A(Row, Stage) - Sum);
           end loop;
	else
           for Row in Stage+1..Final_Index loop
              A(Row, Stage) := Scale_Factor * A(Row, Stage);
           end loop;
        end if;

     end loop; -- Stage

     -- Step 3: Get final row and column.

     Final_Stage := Final_Index;
     Sum := Zero;
     for K in Starting_Index..Final_Stage-1 loop
        Sum := Sum + A(Final_Stage, K) * A(Final_Stage, K);
     end loop;
     Pivot_Info := A(Final_Stage, Final_Stage) - Sum;

     if Pivot_Info < Min_Allowed_Real then
        raise Constraint_Error with "Matrix doesn't seem to be positive definite.";
     end if;

     A(Final_Stage, Final_Stage) := Sqrt (Pivot_Info);
     Diag_Inverse(Final_Stage)   := One / A(Final_Stage, Final_Stage);

     -- Step 4: Have L, now fill in the Upper triangular
     -- matrix with transpose(L).  Recall A_LU = L*U = L*transpose(L).
     -- Read from the Lower triangle.  Row > Col for the Lower Tri:

     for Row in Starting_Index+1 .. Final_Index loop
        for Col in Starting_Index..Row-1 loop
           A(Col, Row) := A(Row, Col);
        end loop;
     end loop;

  end LU_Decompose;

  -----------
  -- Solve --
  -----------
  
  --  Solve for X in the equation A*X = b.  Below you enter the LU
  --  decomposition of A, not A itself.  A_Lu and Diag_Inverse are
  --  the objects returned by LU_decompose.
  
  procedure Solve 
    (X              :    out Col_Vector;
     B              : in     Col_Vector;
     A_LU           : in     Matrix;
     Diag_Inverse   : in     Col_Vector;
     Final_Index    : in     Index := Index'Last;
     Starting_Index : in     Index := Index'First)
  is
     Z : Col_Vector;
     First_Non_Zero_B : Index := Starting_Index;
     Sum              : Real;
  begin
     for Row in Index loop
        X(Row) := Zero;
     end loop;

     --for I in Starting_Index..Final_Index loop
     --   if Abs (B(I)) > 0.0 then
     --      First_Non_Zero_B := I;
     --      exit;
     --   end if;
     --end loop;
     --
     --  In solving for Z in the equation L Z = B, the Z's are
     --  zero up to the 1st non-zero element of B.
     --
     --if First_Non_Zero_B > Starting_Index then
     --   for I in Index range Starting_Index..First_Non_Zero_B-1 loop
     --      Z(I) := 0.0;
     --   end loop;
     --end if;
     
     First_Non_Zero_B := Starting_Index;

     --  The matrix equation is in the form L * tr(L) * X2 = B.
     --  First assume tr(L) * X2 is Z, and
     --  solve for Z in the equation L Z = B.

     Z(First_Non_Zero_B) := B(First_Non_Zero_B) * Diag_Inverse(First_Non_Zero_B);

     if First_Non_Zero_B < Final_Index then
     
        for Row in First_Non_Zero_B+1..Final_Index loop
           Sum := Zero;
           for Col in Index range Starting_Index..Row-1 loop
              Sum := Sum + A_LU(Row, Col) * Z(Col);
           end loop;
           Z(Row) := (B(Row) - Sum) * Diag_Inverse(Row);
        end loop;
        
     end if;

     --  Solve for X2 in the equation tr(L) * X2 = U * X2 = Z.

     X(Final_Index) := Z(Final_Index) * Diag_Inverse(Final_Index);
     
     if Final_Index > Starting_Index then
     
        for Row in reverse Starting_Index .. Final_Index-1 loop
           Sum := Zero;
           for Col in Index range Row+1..Final_Index loop
              Sum := Sum + A_LU(Row, Col) * X(Col);
           end loop;
           X(Row) := (Z(Row) - Sum) * Diag_Inverse(Row);
        end loop;
        
     end if;

  end Solve;
   
end Cholesky_LU;
