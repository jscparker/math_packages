
--------------------------------------------------------------------------
-- package body Jacobi_Eigen, Jacobi iterative eigen-decomposition
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
---------------------------------------------------------------------------

with Hypot;

package body Jacobi_Eigen is

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;
   Half : constant Real := +0.5;

   package Hyp is new Hypot (Real); use Hyp;

   Min_Exp : constant Integer := Real'Machine_Emin;
   Min_Allowed_Real : constant Real := Two ** (Min_Exp - Min_Exp/4);

  ---------------------------------
  -- Get_Jacobi_Rotation_Factors --
  ---------------------------------

  --  underflows OK here. No overflows are OK.

  procedure Get_Jacobi_Rotation_Factors
    (P, Q    : in    Real;
     s       :    out Real;
     tau     :    out Real;
     Del_D   :    out Real)
  is
     t, Alpha, A, B, Denom : Real;          -- t is Q / P
  begin
     s     := Zero;
     tau   := Zero;
     Del_D := Zero;

     if  Abs (P) < Min_Allowed_Real  then
        return; -- use default vals of s, tau, del_d
     end if;

     if  Abs (Q) > Abs (P)  then

        -- Abs (P) > 0 (tested before arrival also), which implies Abs (Q) > 0.

        t := P / Q;

        Alpha := Half * t * (t / (One + Hypotenuse (One, t)));
        A     := Half * t;
        B     := One + Alpha;

     else       -- if Abs (P) >= Abs (Q) then

        t := Q / P;

        Alpha := Abs (t) + t * (t / (One + Hypotenuse (One, t)));
        A     := One;
        B     := One + Alpha;
        if t < Zero then A := -A; end if;

     end if;

   --s       := A / Sqrt(B*B + A*A);                  -- Sine
   --c       := B / Sqrt(B*B + A*A);                  -- Cosine

     Denom   := Hypotenuse (A, B);
     s       := A / Denom;
     tau     := A / (B + Denom);                  --  -(cos - 1) / sin
     Del_D   := (P * A) / B;

  end Get_Jacobi_Rotation_Factors;

  ---------------------
  -- Eigen_Decompose --
  ---------------------

  procedure Eigen_Decompose
    (A                      : in out Matrix;
     Q_tr                   :    out Matrix;
     Eigenvals              :    out Col_Vector;
     No_of_Sweeps_Performed :    out Natural;
     Total_No_of_Rotations  :    out Natural;
     Start_Col              : in     Index   := Index'First;
     Final_Col              : in     Index   := Index'Last;
     Eigenvectors_Desired   : in     Boolean := True)
  is
     D  : Col_Vector renames Eigenvals;
     Z, B  : Col_Vector;

     Max_Allowed_No_of_Sweeps : constant Positive := 256; -- badly scaled A needs lots
     No_of_Preliminary_Sweeps : constant Positive := 14;  -- use 14; don't change.

     subtype Sweep_Range is Positive range 1 .. Max_Allowed_No_of_Sweeps;

     Reciprocal_Epsilon  : constant Real := One / (Real'Epsilon * Two**(-3));
     -- Bst stnd setting for accuracy seems to be about:  Real'Epsilon * Two**(-3).
     -- Matrices with clustered eigvals seem to need the more accurate setting (3).
     -- Usually, Real'Epsilon := 2.0**(-50) for 15 digit Reals.

     Matrix_Size : constant Real := Real (Final_Col) - Real (Start_Col) + One;
     No_of_Off_Diag_Elements : constant Real := Half*Matrix_Size*(Matrix_Size-One);
     Exp : Integer;
     Factor : Real;

     s, g, h, tau : Real;   -- Rutishauser variable names.
     Q, Del_D   : Real;
     Sum, Mean_Off_Diagonal_Element_Size, Threshold : Real;
     Pivot : Real;

  begin

     -- Initialize all out parameters. D renames Eigenvals.
     -- Q_tr starts as Identity; is rotated into set of Eigenvectors of A.

     if Eigenvectors_Desired then
        Q_tr := (others => (others => Zero));
        for j in Index loop
           Q_tr(j, j) := One;
        end loop;
     end if;
     --  Don't fill (potentially) giant array with 0's unless it's needed.
     --  If Eigenvectors_Desired=False, we can free up much memory if
     --  this is never touched.

     Z := (others => Zero);
     B := (others => Zero);
     D := (others => Zero);
     for j in Start_Col .. Final_Col loop  -- assume A not all init
       D(j) := A(j, j);
       B(j) := A(j, j);
     end loop;

     No_of_Sweeps_Performed := 0;
     Total_No_of_Rotations  := 0;

     if Matrix_Size <= One then  return;  end if; -- right answer for Size=1.


     Sweep_Upper_Triangle:
     for Sweep_id in Sweep_Range loop

        No_of_Sweeps_Performed := Sweep_id - Sweep_Range'First;

        Sum := Zero;
        for Row in Start_Col .. Final_Col-1 loop    --sum off-diagonal elements
        for Col in Row+1 .. Final_Col loop
           Sum := Sum + Abs (A(Row, Col));
        end loop;
        end loop;
        Mean_Off_Diagonal_Element_Size := Sum / No_of_Off_Diag_Elements;


        exit Sweep_Upper_Triangle when Mean_Off_Diagonal_Element_Size < Min_Allowed_Real;


        --  Program does Threshold pivoting.
        --
        --  If a Pivot (an off-diagonal matrix element) satisfies
        --  Abs (Pivot) > Threshold, then do a Jacobi rotation to zero it out.
        --
        --  Next calculate size of Threshold:

        if Sweep_id > No_of_Preliminary_Sweeps then
           Threshold := Zero;
        elsif Standard_Threshold_Policy then
           Threshold := One * Mean_Off_Diagonal_Element_Size;
           --  On average, fewer overall rotations done here, at slight
           --  expense of accuracy.
        else
           Exp       := 11 - Sweep_id;
           Factor    := One + (+1.666)**Exp;
           Threshold := Factor * Mean_Off_Diagonal_Element_Size;
           --  The big Threshold here helps with badly scaled matrices.
           --  May improve accuracy a bit if scaling is bad. Policy
           --  here is closer to that of the original Jacobi, which
           --  always rotates away the largest Pivots 1st.
        end if;

        Pivots_Row_id: for Pivot_Row in Start_Col .. Final_Col-1 loop
        sum := 0.0;
        Pivots_Col_id: for Pivot_Col in Pivot_Row+1 .. Final_Col loop

           Pivot := A(Pivot_Row, Pivot_Col);

           --  Have to zero-out sufficiently small A(Pivot_Col, Pivot_Row) to get convergence, 
           --  ie, to get  Mean_Off_Diagonal_Element_Size -> 0.0. The test is:
           --
           --  A(Pivot_Col, Pivot_Row) / Epsilon <= Abs D(Pivot_Col) and 
           --  A(Pivot_Col, Pivot_Row) / Epsilon <= Abs D(Pivot_Row).

           if (Sweep_id > No_of_Preliminary_Sweeps) and then
              (Abs (Pivot) * Reciprocal_Epsilon <= Abs D(Pivot_Row)) and then
              (Abs (Pivot) * Reciprocal_Epsilon <= Abs D(Pivot_Col))     then

              A(Pivot_Row, Pivot_Col) := Zero;

           elsif Abs (Pivot) > Threshold then

              Q := Half * (D(Pivot_Col) - D(Pivot_Row));

              Get_Jacobi_Rotation_Factors (Pivot, Q, s, tau, Del_D);

              D(Pivot_Row) := D(Pivot_Row) - Del_D;  -- Locally D is only used for threshold test.
              D(Pivot_Col) := D(Pivot_Col) + Del_D;
              Z(Pivot_Row) := Z(Pivot_Row) - Del_D;  -- Z is reinitialized to 0 each sweep, so
              Z(Pivot_Col) := Z(Pivot_Col) + Del_D;  -- it sums the small d's 1st. Helps a tad.

              A(Pivot_Row, Pivot_Col) := Zero;

              if Pivot_Row > Start_Col then
              for j in Start_Col .. Pivot_Row-1 loop
                 g := A(j, Pivot_Row);
                 h := A(j, Pivot_Col);
                 A(j, Pivot_Row) := g-s*(h+g*tau);
                 A(j, Pivot_Col) := h+s*(g-h*tau);
              end loop;
              end if;

              for j in Pivot_Row+1 .. Pivot_Col-1 loop
                 g := A(Pivot_Row, j);
                 h := A(j, Pivot_Col);
                 A(Pivot_Row, j) := g-s*(h+g*tau);
                 A(j, Pivot_Col) := h+s*(g-h*tau);
              end loop;

              if Pivot_Col < Final_Col then
              for j in Pivot_Col+1 .. Final_Col loop
                 g := A(Pivot_Row, j);
                 h := A(Pivot_Col, j);
                 A(Pivot_Row, j) := g-s*(h+g*tau);
                 A(Pivot_Col, j) := h+s*(g-h*tau);
              end loop;
              end if;

              if Eigenvectors_Desired then
              for j in Start_Col .. Final_Col loop
                 g := Q_tr(Pivot_Row, j);
                 h := Q_tr(Pivot_Col, j);
                 Q_tr(Pivot_Row, j) := g-s*(h+g*tau);
                 Q_tr(Pivot_Col, j) := h+s*(g-h*tau);
              end loop;
              end if;

              Total_No_of_Rotations := Total_No_of_Rotations + 1;

           end if; -- if (Sweep_id > No_of_Preliminary_Sweeps)

        end loop Pivots_Col_id;
        end loop Pivots_Row_id;

        for j in Start_Col .. Final_Col loop  -- assume A not all initialized
           B(j) := B(j) + Z(j);
           D(j) := B(j);
           Z(j) := Zero;
        end loop;

     end loop Sweep_Upper_Triangle;

  end Eigen_Decompose;

  ---------------
  -- Sort_Eigs --
  ---------------

  procedure Sort_Eigs
    (Eigenvals         : in out Col_Vector;
     Q_tr              : in out Matrix;             -- rows  are the eigvectors
     Start_Col         : in     Index   := Index'First;
     Final_Col         : in     Index   := Index'Last;
     Sort_Eigvecs_Also : in     Boolean := True)
  is
     Max_Eig, tmp : Real;
     Max_id  : Index;
  begin
     if Start_Col < Final_Col then
     for i in Start_Col .. Final_Col-1 loop

        Max_Eig := Eigenvals(i);   Max_id := i;

        for j in i+1 .. Final_Col loop
           if Eigenvals(j) > Max_Eig then
              Max_Eig := Eigenvals(j);   Max_id := j;
           end if;
        end loop;

        tmp               := Eigenvals(i);
        Eigenvals(i)      := Max_Eig;
        Eigenvals(Max_id) := tmp;

        -- swap rows of Q_tr:

        if Sort_Eigvecs_Also then
        for k in Start_Col .. Final_Col loop
           tmp             := Q_tr(i, k);
           Q_tr(i, k)      := Q_tr(Max_id, k);
           Q_tr(Max_id, k) := tmp;
        end loop;
        end if;

     end loop;
     end if;

  end Sort_Eigs;

end Jacobi_Eigen;

