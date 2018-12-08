
---------------------------------------------------------------------------
-- package body QR_Symmetric_Eigen, QR based eigen-decomposition
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

with Ada.numerics.Generic_Elementary_Functions;
with Tridiagonal;
with Hypot;

package body QR_Symmetric_Eigen is

  package Maths is new Ada.numerics.Generic_Elementary_Functions(Real); use Maths;

  package Hypo is new Hypot (Real);

  package Tridi is new Tridiagonal (Real, Index, Matrix);

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  --
  -- Matrix
  --
  --   |a  b|
  --   |b  c|
  --
  -- has eigenvals  x  satisfying:
  --
  --   x^2  - (a + c) * x  +  a*c - b*b  =  0
  --   x^2  - (a + c) * x  +  a_0        =  0
  --
  -- Solution (x) is always real:
  --
  --   x = (a + c)/2 +/- Sqrt (((a + c)/2)^2 - a_0)
  --   x = (a + c)/2 +/- Sqrt (((a - c)/2)^2 + b*b)
  --
  --   Let y = x - c, and Beta = (a - c)/2. Then
  --
  --   y = (a - c)/2 +/- Sqrt (((a - c)/2)^2 + b*b)
  --   y =  Beta     +/- Sqrt (Beta^2 + b*b)
  --
  -- General Soln:
  --
  --   Eig_Val - c = y = Beta * (1 +/- Sgn(Beta) * Sqrt (Beta^2 + b*b))
  --
  -- Use
  --
  --   Sqrt (Beta^2 + b*b) - 1 = b^2 / (Beta^2 + |Beta| * Sqrt (Beta^2 + b*b))
  --
  -- and
  --
  --   1 + Sqrt (Beta^2 + b*b) = 2 + (Sqrt (Beta^2 + b*b) - 1)
  --
  -- to get
  --
  --   Eig_1 = a + Sgn(Beta) * b^2 / (|Beta| + Sqrt (Beta^2 + b*b))
  --   Eig_2 = c - Sgn(Beta) * b^2 / (|Beta| + Sqrt (Beta^2 + b*b))
  --
  --
  procedure Quadratic_Eigs
    (a, b, c      : in Real;
     Eig_Val_1    : out Real;
     Eig_Val_2    : out Real)
  is
     Beta : constant Real := 0.5 * (a - c);
  begin
  
     if  b = Zero  then
        Eig_Val_1 := a;
        Eig_Val_2 := c;
        return;
     end if;

     if  Beta = Zero  then   -- a = c
        Eig_Val_1 := a + Abs b;
        Eig_Val_2 := c - Abs b;
        return;
     end if;

     -- Hypot = Sqrt (Beta^2 + b^2);

     declare
        Hypot : constant Real := Hypo.Hypotenuse (b, Beta);
        Shift : Real := b * (b / (Hypot + Abs Beta));
     begin
        Shift     := Real'Copy_Sign (Shift, Beta);
        Eig_Val_1 := a + Shift;
        Eig_Val_2 := c - Shift;
     end;
 
     --     --   test using:  x^2  - (a + c) * x  +  a_0  =  0
     --  
     --     declare
     --        x, tst : Real;
     --        Min_R : Real := Min_Allowed_Real;
     --        a_0   : constant Real := a*c - b*b;
     --     begin
     --        x := Eig_Val_1;
     --        tst := ((x - a)*(x - c) - b*b) / (x*x + b*b + Abs a_0 + Min_R);
     --        if Abs tst > Real'Epsilon * 16.0 then 
     --           new_line; put(Real'Image(tst)); new_line;
     --           null;
     --        end if;
     --        x := Eig_Val_2;
     --        tst := ((x - a)*(x - c) - b*b) / (x*x + b*b + Abs a_0 + Min_R);
     --        if Abs tst > Real'Epsilon * 16.0 then 
     --           new_line; put(Real'Image(tst)); new_line;
     --           null;
     --        end if;
     --     end;
  
  end Quadratic_Eigs;

  procedure Swap (A, B : in out Real) is 
     tmp : constant Real := A;
  begin
     A := B; B := tmp;
  end Swap;

  ---------------
  -- Sort_Eigs --
  ---------------

  procedure Sort_Eigs
    (Eigenvals         : in out Col_Vector;
     Q                 : in out Matrix;             -- Columns are the eigvectors
     Start_Col         : in     Index   := Index'First;
     Final_Col         : in     Index   := Index'Last;
     Sort_Eigvecs_Also : in     Boolean := False)
  is
     Max_Eig : Real;
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

        -- swap i with Max_id:

        Swap (Eigenvals(i), Eigenvals(Max_id));

        if Sort_Eigvecs_Also then
        for r in Start_Col .. Final_Col loop
           Swap (Q(r, i), Q(r, Max_id));
        end loop;
        end if;

     end loop;
     end if;

  end Sort_Eigs;

  ----------
  -- Norm --
  ----------

  function Norm 
    (Q            : in Col_Vector;
     Starting_Col : in Index := Index'First;
     Final_Col    : in Index := Index'Last)
     return Real
  is
     Max : Real := Zero;
     Sum, tst : Real := Zero;
     Scale_Factor  : Real := One;
     Scaled : Boolean := False;
     Scale_Exp : Integer;
     Max_Exp               : constant Integer := Real'Machine_Emax - 8;
     Sqrt_Max_Allowed_Real : constant Real    := Two**(Max_Exp / 2);
     Min_Exp               : constant Integer := Real'Machine_Emin + 8;
     Sqrt_Min_Allowed_Real : constant Real    := Two**(Min_Exp / 2);
  begin
     for i in Starting_Col .. Final_Col  loop
        tst := Abs (Q(i));
        if tst > Max then Max := tst; end if;
     end loop;

     if Max > Sqrt_Max_Allowed_Real then
        Scale_Exp    := Integer'Min (Real'Exponent (Max), Real'Machine_Emax-4);
        Scale_Factor := Two**(-Scale_Exp);
        Scaled       := True;
     end if;

     if Max < Sqrt_Min_Allowed_Real then
        Scale_Exp    := Integer'Max (Real'Exponent (Max), Real'Machine_Emin+4);
        Scale_Factor := Two**(-Scale_Exp);
        Scaled       := True;
     end if;

     if Scaled then
        for k in Starting_Col .. Final_Col loop
           Sum := Sum + (Scale_Factor * Q(k))**2;
        end loop;
     else
        for k in Starting_Col .. Final_Col loop
           Sum := Sum + Q(k)**2;
        end loop;
     end if;

     if Scaled then
        return Sqrt (Abs Sum) * Two**Scale_Exp;
     else
        return Sqrt (Abs Sum);
     end if;

  end Norm;

  ---------------------
  -- Eigen_Decompose --
  ---------------------

  procedure Eigen_Decompose
    (A                      : in out Matrix;
     Q                      :    out Matrix;
     Eigenvals              :    out Col_Vector;
     Start_Col              : in     Index   := Index'First;
     Final_Col              : in     Index   := Index'Last;
     Eigenvectors_Desired   : in     Boolean := False)
  is
     Shift, e1, e2 : Real;
     N : constant Integer := Integer (Final_Col) - Integer (Start_Col) + 1;
     Min_Exp : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real := 2.0 ** (Min_Exp - Min_Exp/32);
  begin
  
     -- Q starts as Identity; is rotated into set of Eigenvectors of A.

     Q := (others => (others => Zero));
     for j in Index loop
        Q(j,j) := One;
     end loop;

     Eigenvals := (others => Zero);
   
     --  A_tridiag = Q_tr * A_true * Q.
     --
     --  A_true    =  Q * A_tridiag * Q_tr
  
     Tridi.Tridiagonalize
       (A                => A,        -- A = A_true is replaced with:  A_tridiag
        Q                => Q,        -- A_true  =  Q * A_tridiag * Q_tr
        Starting_Col     => Start_Col,
        Final_Col        => Final_Col,
        Initial_Q        => Q,
        Q_Matrix_Desired => Eigenvectors_Desired);
  
     -- A is now in Upper_Hessenberg form, thanks to Q:
     --
     --    A_true = Q A Q' = Q A_tridiag Q'
     --
     -- Procedure Upper_Hessenberg initialized "out" parameter Q as the new Z_r.
     --
     for i in reverse Start_Col+1 .. Final_Col loop

        Iterate: for iter_id in 0 .. 30 + N / 4 loop  -- 15 iters is very rare

        if Abs A(i,i-1) < Min_Allowed_Real then
           A(i,i-1) := Zero; A(i-1,i) := Zero;
        end if;

        --if Abs A(i,i-1) = Zero then
        --   text_io.put(integer'image(iter_id)); exit;
        --end if;
        exit Iterate when Abs A(i,i-1) = Zero;

        Quadratic_Eigs
          (a          => A(i-1,i-1), 
           b          => A(i,i-1), 
           c          => A(i,i),
           Eig_Val_1  => e1,
           Eig_Val_2  => e2);

        if Abs(e2-A(i,i))<Abs(e1-A(i,i)) then Shift:=e2; else Shift:=e1; end if;

        Tridi.Lower_Diagonal_QR_Iteration
          (A                => A,
           Q                => Q,
           Shift            => Shift,
         --Final_Shift_Col  => Final_Col,
           Final_Shift_Col  => i,  -- short cut, hardly matters
           Starting_Col     => Start_Col,
           Final_Col        => Final_Col,
           Q_Matrix_Desired => Eigenvectors_Desired);

        end loop iterate;
     end loop;

     -- Eigenvectors of A are returned as the Columns of matrix Q.
     --
     -- so   Q_tr * A * Q  =  D  =  Diagonal_Eigs
     --
     -- so   Q * D * Q_tr  =  A  =  Original Matrix

     for i in Index range Start_Col .. Final_Col loop 
        Eigenvals(i) := A(i,i); 
     end loop;

  end Eigen_Decompose;

end QR_Symmetric_Eigen;

