
---------------------------------------------------------------------------
-- package body Givens_Rotation
-- Copyright (C) 2018 Jonathan S. Parker.
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

package body Givens_Rotation is

  package Hypotenuse is new Hypot (Real); use Hypotenuse;

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  -----------------------------------
  -- Get_Rotation_That_Zeros_Out_L --
  -----------------------------------

  -- P = Pivot, L = Low.
  --
  -- cos = P/r,   sin = L/r,   Hypot = r = sqrt(P*P + L*L)
  --
  -- clockwise rotation: notice all rotations are centered on the diagonal
  --
  --   1  0  0  0  0       0       0
  --   0  c  s  0  0       P       r
  --   0 -s  c  0  0   x   L   =   0
  --   0  0  0  1  0       0       0
  --   0  0  0  0  1       0       0
  --
  --
  -- if |L| >= |P| then   t = P/L  <= 1
  -- if |P| >  |L| then   t = L/P  <  1
  --
  --
  --  let t = smaller / larger
  --
  --  let u = 1/sqrt(1+t*t) and w = sqrt(t*t/(1+t*t)) with
  --
  --  use (1 - 1/sqrt(1+t*t)) * (1 + 1/sqrt(1+t*t)) = t*t / (1 + t*t)
  --
  --  1/sqrt(1+t*t) - 1 = - t*t/(sqrt(1+t*t) + 1+t*t)
  --                    = -(t/(sqrt(1+t*t))*(t/(1 + Sqrt(1+t*t)))
  --                    = - Abs (w) * Abs (t)/(1+ sqrt(1+t*t))
  --                    =   u_lo
  --
  --  u_hi = 1   =>   u_lo + u_hi  =  1/sqrt(1+t*t) = u = Abs (cs) if a<b
  --
  --  hypot = |L| * sqrt(1+t*t) = |L| * (1  + t*t/(1+sqrt(1+t*t)) = hi + lo if t<1
  --
  --   P : Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
  --   L : Real := A(Low_Row, Pivot_Col);

  procedure Get_Rotation_That_Zeros_Out_Low
    (Pivot, Low      : in     Real;
     sn, cs          :    out Real;
     cs_minus_1      :    out Real;
     sn_minus_1      :    out Real;
     hypot           :    out Real;
     P_bigger_than_L :    out Boolean;
     Skip_Rotation   :    out Boolean)
  is
     Emax : constant Integer := Real'Machine_Emax;
     Emin : constant Integer := Real'Machine_Emin;
     P, L, Abs_P, Abs_L : Real;
     min_arg_over_hypot, max_arg_over_hypot_minus_1 : Real;
  begin

     if (not Low'Valid) or (not Pivot'Valid) then
        raise Constraint_Error with "Invalid input data in Get_Rotation...";
     end if;

     P := Pivot; L := Low;

     Abs_P := Abs (P); Abs_L := Abs (L);

     -- default: no rotation is performed.
     sn              := Zero;
     cs              := One;
     sn_minus_1      := -One;
     cs_minus_1      := Zero;
     hypot           := Abs_P;
     P_bigger_than_L := True;
     Skip_Rotation   := True;

     if Abs_L < Two**Emin then 
        Skip_Rotation := True; -- use defaults
        return;
     end if;

     if Abs_P < Two**Emin then
        if Abs_L > Two**Emin then 
           sn              := One;
           cs              := Zero;
           sn_minus_1      := Zero;
           cs_minus_1      := -One;
           hypot           := Abs_L;
           P_bigger_than_L := False;
           Skip_Rotation   := False;
           return;
        else
           Skip_Rotation := True; -- use defaults
           return;
        end if;
     end if;

     -- if the Low val is too small compared to the pivot, then
     -- it contributes nothing to pivot after rotation ... but the row
     -- still might contribute elsewhere.

     if Abs_L < Two**(Emax - Emax/10 - Real'Machine_Radix - 1) and then 
        Abs_L * Two**(Emax/10 + Real'Machine_Radix) < Abs_P
     then 
        Skip_Rotation := True;
        return;  -- use defaults, Skip_Rotation.
     end if; -- essential optimisation.

     if Abs_P > Abs_L then   -- |s| <= |c|

          Skip_Rotation   := False;
          P_bigger_than_L := True;

          -- Want P>0 always, (for efficiency, since we use cos = P/R = 1 + (cos-1)),
          -- so if P<0 then flip signs of both P and L
          -- to ensure zero'd out element.
          if P < Zero then  P := -P; L := -L; end if;

          Get_Hypotenuse (P, L, hypot, min_arg_over_hypot, max_arg_over_hypot_minus_1);

        --cs                 := P / hypot;  -- normally unused; set to default

          sn                 :=  Real'Copy_Sign (min_arg_over_hypot, L);
        --sn                 :=  L / hypot;

          cs_minus_1         :=  max_arg_over_hypot_minus_1;
        --cs_minus_1         := -Abs (sn) * Abs_L / hypot_plus_max_arg;
        --cs_minus_1_over_sn := -L / hypot_plus_max_arg;

     else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1: NOTICE <= !

          Skip_Rotation   := False;
          P_bigger_than_L := False;

          -- Want L>0 always. If L<0 then flip signs of both P and L
          -- to ensure zero'd out element.
          if L < Zero then  P := -P; L := -L; end if;

          Get_Hypotenuse
            (P, L, hypot, min_arg_over_hypot, max_arg_over_hypot_minus_1);

        --sn                 :=  L / hypot;        -- set to default; unused

          cs                 :=  Real'Copy_Sign (min_arg_over_hypot, P);
        --cs                 :=  P / hypot;

          sn_minus_1         :=  max_arg_over_hypot_minus_1;
        --sn_minus_1         := -Abs (cs) * Abs_P / hypot_plus_max_arg;
        --sn_minus_1_over_cs := -P / hypot_plus_max_arg;

     end if;

  end Get_Rotation_That_Zeros_Out_Low;

end Givens_Rotation;

