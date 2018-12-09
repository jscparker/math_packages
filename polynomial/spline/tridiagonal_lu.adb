
-----------------------------------------------------------------------
-- package body Tridiagonal_LU, LU decomposition for Tridiagonal matrices.
-- Copyright (C) 2018 Jonathan S. Parker
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
---------------------------------------------------------------------------

package body Tridiagonal_LU is

  -- The upper matrix is U, the lower L.
  -- Assume that the diagonal elements of U are 1.0. So
  -- the diagonal elements of L but not U appear on the
  -- the diagonal of the output matrix A.

  procedure LU_Decompose
    (A            : in out Matrix;
     Index_Start  : in     Index := Index'First;
     Index_Finish : in     Index := Index'Last)
  is
     Stage : Index;
  begin
     if Index_Finish < Index_Start then return; end if;

     -- Stage 1:
     Stage := Index_Start;

     if Abs (A(0)(Stage)) < Epsilon then
        if Set_Zero_Valued_Pivots_To_Epsilon then
          A(0)(Stage) := Epsilon;
        else
          raise matrix_is_singular;
        end if;
     end if;

     A(1)(Stage) := A(1)(Stage) / A(0)(Stage);

     for Stage in Index_Start+1 .. Index_Finish-1 loop

        -- At each stage we calculate row "stage" of the Upper matrix U
        -- and Column "Stage" of the Lower matrix L.
        -- The matrix A is overwritten with these, because the elements
        -- of A in those places are never needed in future stages.
        -- However, the elements of U and L ARE needed in those places,
        -- so to get those elements we will be accessing A (which stores them).
        --
        -- Get row "stage" of U and column "stage" of L. Notice these formulas
        -- update only (Stage+1..Last) elements of the respective row
        -- and column, and depend on only (1..Stage) elements
        -- of U and L, which were calculated previously, and stored in A.

        A(0)(Stage) := A(0)(Stage) - A(-1)(Stage)*A(1)(Stage-1);

        if Abs (A(0)(Stage)) < Epsilon then
           if Set_Zero_Valued_Pivots_To_Epsilon then
              A(0)(Stage) := Epsilon;
           else
              raise matrix_is_singular;
           end if;
        end if;

        A(1)(Stage) := A(1)(Stage) / A(0)(Stage);

     end loop;

      -- Step 3: final row and column.

     if Index_Finish > Index_Start then
        Stage := Index_Finish;
        A(0)(Stage) := A(0)(Stage) - A(-1)(Stage)*A(1)(Stage-1);
     end if;

  end LU_Decompose;


  procedure Solve
    (X            :    out Column;
     A            : in     Matrix;
     B            : in     Column;
     Index_Start  : in     Index := Index'First;
     Index_Finish : in     Index := Index'Last)
  is
     Z2               : Column;
     FirstNonZeroB    : Index := Index_Start;
     Xtemp, Xprevious : Real;
  begin
     if Index_Finish < Index_Start then return; end if;

     -- An optimization to make matrix inversion efficient: if B
     -- is a unit vector (all zeros except for a 1.0) then find
     -- 1st non-zero element of B:

     for i in Index range Index_Start..Index_Finish loop
        if Abs (B(i)) > 0.0 then
           FirstNonZeroB := I;
           exit;
        end if;
     end loop;

     -- When solving for Z2 in the equation L Z2 = B, the Z2's will
     -- all be zero up to the 1st non-zero element of B.

     if FirstNonZeroB > Index_Start then
        for i in Index range Index_Start..FirstNonZeroB-1 loop
           Z2(i) := 0.0;
        end loop;
     end if;

     -- Matrix equation is in the form L*U*X2 = B.
     -- Assume U*X2 is Z2, and solve for Z2 in the equation L Z2 = B.
     -- The matrix A contains along its diagonal the
     -- diagonal elements of L - this time remember to divide:

     Z2(FirstNonZeroB) := B(FirstNonZeroB) / A(0)(FirstNonZeroB);

     if FirstNonZeroB < Index_Finish then
        for i in Index range FirstNonZeroB + 1 .. Index_Finish loop
           Z2(i) := (B(i) - A(-1)(i) * Z2(I-1)) / A(0)(i);
        end loop;
     end if;

     -- Solve for X in the equation U X = Z2.
     -- By assumption U = 1.0 on diagonal.

     X(Index_Finish) := Z2(Index_Finish);
     Xprevious       := Z2(Index_Finish);

     for i in reverse Index range Index_Start..Index_Finish-1 loop
        Xtemp     := (Z2(i) - A(1)(i) * Xprevious);
        X(i)      := Xtemp;
        Xprevious := Xtemp;
     end loop;

  end Solve;
end Tridiagonal_LU;
