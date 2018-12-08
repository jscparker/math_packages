
----------------------------------------------------------------------
-- package body Peters_Eigen, eigendecomposition
-- Copyright (C) 2008-2018 Jonathan S. Parker.
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
with Hessenberg;
with Givens_QR_Method;
with Hypot;

package body Peters_Eigen is

  package Givens_Hess is new Hessenberg (Real, Index, Matrix);
  package Initial_QR  is new Givens_QR_Method (Real, Index, Matrix);
  package Hypo        is new Hypot (Real);

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;
  Half : constant Real := +0.5;

  Max_Exp               : constant Integer := Real'Machine_Emax - Real'Machine_Emax/16;
  Sqrt_Max_Allowed_Real : constant Real    := Two**(Max_Exp / 2);

  Min_Exp               : constant Integer := Real'Machine_Emin - Real'Machine_Emin/16;
  Sqrt_Min_Allowed_Real : constant Real    := Two**(Min_Exp / 2);

  package Math is new Ada.numerics.Generic_Elementary_Functions(Real);
  use Math;

  ----------
  -- Div2 --
  ----------

  -- Complex division, C := A / B
  --
  -- If B = 0, allow the usual real arithmetic to raise exceptions.
  --
  -- Algorithm by Robert L. Smith. (Supposed to be better at avoiding
  -- unnecessary divisions by 0 (due to underflow) than Eispack version.)
  --
  procedure Div2
    (A_r, A_i : in     Real;
     B_r, B_i : in     Real;
     C_r, C_i :    out Real)
  is
     e, f : Real;
  begin
     if Abs (B_i) < Abs (B_r) then
        e   :=  B_i / B_r;
        f   :=  B_r + B_i*e;
        C_r := (A_r + A_i*e) / f;
        C_i := (A_i - A_r*e) / f;
     else
        e   :=  B_r / B_i;
        f   :=  B_i + B_r*e;
        C_r := (A_i + A_r*e) / f;
        C_i :=(-A_r + A_i*e) / f;
     end if;
  end Div2;

  procedure Div
    (A_re, A_im : in     Real;
     B_re, B_im : in     Real;
     C_r, C_i   :    out Real)
  is
     A_max : constant Real := Real'Max (Abs A_re, Abs A_im);
     B_max : constant Real := Real'Max (Abs B_re, Abs B_im);
     A_Scale, B_Scale, Unscale, B_Squ : Real := One;
     A_Scale_exp, B_Scale_exp : Integer := 0;
     A_r, A_i, B_r, B_i : Real;
  begin
     C_r := Zero;
     C_i := Zero;
     if A_max = Zero then -- use default C_r, C_i
        return;
     end if;

     if B_max = Zero then
        raise Constraint_Error with "Division by 0 in Div";
     end if;

     A_r := A_re; A_i := A_im;
     if A_max > Two**(Real'Machine_Emax/2 - 2) or else A_max < Two**(Real'Machine_Emin/2 + 2) then
        A_Scale_exp := -Real'Exponent (A_max);
        A_Scale     := Two**A_Scale_exp;
        A_r := A_re * A_Scale;
        A_i := A_im * A_Scale;
     end if;

     B_r := B_re; B_i := B_im;
     if B_max > Two**(Real'Machine_Emax/2 - 2) or else B_max < Two**(Real'Machine_Emin/2 + 2) then
        B_Scale_exp := -Real'Exponent (B_max);
        B_Scale     := Two**B_Scale_exp;
        B_r := B_re * B_Scale;
        B_i := B_im * B_Scale;
     end if;

     Unscale := Two**(B_Scale_exp - A_Scale_exp);
     -- over/under flow if A/B out of range.

     B_Squ := B_r**2 + B_i**2;
     C_r := ((A_r*B_r + A_i*B_i) / B_Squ) * Unscale;
     C_i := ((A_i*B_r - A_r*B_i) / B_Squ) * Unscale;

  end Div;

  --------------------
  -- Balance_Matrix --
  --------------------

  -- procedure Balance_Matrix returns Vector D representing a
  -- diagonal matrix, and Matrix B,
  --
  -- B :=  D^(-1) *  A * D
  --
  -- If we subsequently find the eigendecomposition of B
  --
  --    B*v = lambda*v.
  --
  -- Then
  --
  --     D^(-1)*A*D*v = lambda*v
  --
  -- implies that
  --
  --     A*D*v = lambda*D*v,
  --
  -- so that the eigenvectors of A are D*v, (D times the eigenvectors of B).
  -- The eigenvalues A and B are unchanged by the transformation.
  --
  -- The eigenvectors v of matrix B are orthonormal, but not the
  -- eigenvectors D*v of matrix A.
  --
  -- A is written over with B.

  procedure Balance_Matrix 
    (A              : in out Matrix; 
     D              : in out Col_Vector;
     Starting_Col   : in     Index;
     Final_Col      : in     Index;
     Balance_Policy : in     Balance_Code := Partial)
  is
     Radix                    : constant Real := Two ** 1;  -- Two**1 most accurate.
     Radix_Squared            : constant Real := Radix ** 2;
     Reciprocal_Radix         : constant Real := 1.0 / Radix;
     Reciprocal_Radix_Squared : constant Real := Reciprocal_Radix ** 2;

     Max_Allowed_Iterations : Integer := 16;
  
     Converged : Boolean := True;
  
     g, f, s : Real;
     Row_Norm, Col_Norm : Real;
  
     Min_Exp          : constant Integer := Real'Machine_Emin + 16;
     Min_Allowed_Real : constant Real    := Two ** Min_Exp;
  
     Max_Test : constant Real := Two**(Real'Machine_Emax / 2);
     Min_Test : constant Real := Two**(Real'Machine_Emin / 2);
 
     N : constant Integer := 1 + Integer (Real(Final_Col) - Real(Starting_Col));
 
  begin

     case Balance_Policy is
     when Full =>
        Max_Allowed_Iterations := 16;
     when Partial =>
        Max_Allowed_Iterations := 2 + Integer (N / 500);
     when Disabled =>
        Max_Allowed_Iterations := 0;
     end case;
  
     --  Balance a square block whose 
     --        upper-left  corner is (Starting_Col, Starting_Col).
     --        Lower-right corner is (Final_Col, Final_Col).
  
     D := (others => One);
  
     Until_Converged:
     for Interation_id in 1 .. Max_Allowed_Iterations loop
  
        Converged := True;
  
        for i in Starting_Col .. Final_Col loop
     
           Row_Norm := 0.0;
           Col_Norm := 0.0;
           for j in Starting_Col .. Final_Col loop
              if j /= i then
                 Col_Norm := Col_Norm + Abs (A(j, i));
                 Row_Norm := Row_Norm + Abs (A(i, j));
              end if;  
           end loop;
     
           -- D(i) is already 1.0.  If both norms are > 0.0
           -- then calculate f for D(i) and g = 1/f for D(i)^(-1):

           if (Col_Norm > Min_Allowed_Real) and (Row_Norm > Min_Allowed_Real) then 

              f := 1.0;
              g := Row_Norm * Reciprocal_Radix;
              s := Col_Norm + Row_Norm;

              while Col_Norm < g loop
                 exit when f > Max_Test or else Col_Norm > Max_Test;
                 f        := f * Radix;
                 Col_Norm := Col_Norm * Radix_Squared;
              end loop;
     
              g := Row_Norm * Radix;
     
              while Col_Norm > g loop
                 exit when f < Min_Test or else Col_Norm < Min_Test;
                 f        := f * Reciprocal_Radix;
                 Col_Norm := Col_Norm * Reciprocal_Radix_Squared;
              end loop;
    
              -- If we flunk the following test for all row/cols i, then
              -- Converged remains True, and we exit the loop and return.
   
              if ((Row_Norm + Col_Norm) / f  <  0.93 * s) then

                 -- .75 with 3 iters is like .95 with 2 iters

                 Converged := False;
                 D(i) := D(i) * f;            -- D(i) initially 1.0
                 g    := One / f;
     
                 --  (multiplying by D on the right): Scale i_th col by f.
     
                 for k in Starting_Col .. A'Last(1) loop  -- k .. N
                    A(k,i) := f * A(k,i);
                 end loop;
     
                 --  (multiplying by D^(-1) on the left)): Scale i_th row by g = 1/f
     
                 for k in A'First(2) .. Final_Col loop    -- 1 .. l
                    A(i,k) := g * A(i,k);
                 end loop;
              end if;
  
           end if; -- if both norms are non-zero
     
        end loop; -- in i_th row/col
     
        exit when Converged;

     end loop Until_Converged;

  end Balance_Matrix;
 
  ------------------
  -- Balance_Undo --
  ------------------
  --
  -- The eigenvectors of A are D*z, (multiply diagonal matrix
  -- D on the right times the eigenvectors Z of the balanced matrix B).
  --
  procedure Balance_Undo 
    (D            : in     Col_Vector;
     Z            : in out Matrix;
     Starting_Col : in     Index;
     Final_Col    : in     Index)
  is
  begin
   
     for r in Starting_Col .. Final_Col loop
        for c in Starting_Col .. Final_Col loop
           Z(r,c) := Z(r,c) * D(r);
        end loop;
     end loop;
   
  end Balance_Undo;

  ------------------------
  -- Rotation_Factors_2 --
  ------------------------

  --  sn = a / sqrt(a**2 + b**2)
  --  cs = b / sqrt(a**2 + b**2)

  procedure Rotation_Factors_2
    (a, b   : in     Real;
     sn_lo, sn_hi, cs_lo, cs_hi : out Real;
     hypot  :    out Real)
  is
     Abs_a : constant Real := Abs a;
     Abs_b : constant Real := Abs b;
     min_over_hypot, max_over_hypot_minus_1 : Real;
     Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin + 16);
  begin
     --  default:
     cs_hi := One;
     cs_lo := Zero;
     sn_hi := Zero;
     sn_lo := Zero;
     hypot := One;

     if Abs_a >= Abs_b then
        if Abs_a > Min_Allowed_Real then
           Hypo.Get_Hypotenuse (a, b, hypot, min_over_hypot, max_over_hypot_minus_1);
           cs_hi := Zero;
           cs_lo := min_over_hypot;
           sn_hi := One;
           sn_lo := max_over_hypot_minus_1;
           if a < Zero then sn_hi := -sn_hi; sn_lo := -sn_lo; end if;
           if b < Zero then cs_lo := -cs_lo; end if;
        end if;
     else
        if Abs_b > Min_Allowed_Real then
           Hypo.Get_Hypotenuse (a, b, hypot, min_over_hypot, max_over_hypot_minus_1);
           sn_hi := Zero;
           sn_lo := min_over_hypot;
           cs_hi := One;
           cs_lo := max_over_hypot_minus_1;
           if a < Zero then sn_lo := -sn_lo; end if;
           if b < Zero then cs_hi := -cs_hi; cs_lo := -cs_lo; end if;
        end if;
     end if;

  end Rotation_Factors_2;

  ---------------------------------------
  -- Arg_1_is_Negligible_Respect_Arg_2 --
  ---------------------------------------

  -- IF Abs_x is negligible in comparison to Abs_y  then  return True.

  function Arg_1_is_Negligible_Respect_Arg_2 (x, y : Real) return Boolean is
     Abs_x      : constant Real := Abs x;
     Abs_y      : constant Real := Abs y;
     Min_Allowed_Real : constant Real := Two ** Real'Machine_Emin;

     Added_Range        : constant         := 5; -- 2 to 8. 6,7 std. 7 helps with clustered eigs.
     Effective_Mantissa : constant Integer := Real'Machine_Mantissa + Added_Range;
     Eps_Factor         : constant Real    := Two**(-Effective_Mantissa);
     pragma Assert (Added_Range > 1);
     --  Stnd setting for Eps_Factor is:  2**(-Real'Machine_Mantissa-2)
     --  Usually Real'Machine_Mantissa = 53 if Real is 15 digits.
  begin
     if Abs_x < Min_Allowed_Real then -- eg, Abs_x = 0.0
        return True;
     elsif Abs_x <= Abs_y*Eps_Factor then -- Abs_x is negligible in comparison to Abs_y
        return True;
     else
        return False;
     end if;
  end Arg_1_is_Negligible_Respect_Arg_2;

  ------------------------
  -- Unpack_Eig_Vectors --
  ------------------------

  procedure Unpack_Eig_Vectors
    (W_i          : in     Col_Vector;
     Z            : in out Matrix;     -- will exit with Z_r, the real part of Z.
     Z_i          :    out Matrix;     -- will exit with imaginary part.
     Starting_Col : in     Index;
     Final_Col    : in     Index)
  is
  begin

     Z_i := (others => (others => Zero));

     -- Im (W) > 0  at column c means that columns c and c+1 of Z are the
     -- Re and Im parts of the Eigenvector, with eig_val = (W_r(c), Abs(W_i(c))).

     for c in Starting_Col .. Final_Col-1 loop
        if W_i(c) > Zero then
           for j in Starting_Col .. Final_Col loop
              Z_i(j,c)   :=  Z(j,c+1);
              Z_i(j,c+1) := -Z(j,c+1);
           end loop;
           for j in Starting_Col .. Final_Col loop
              Z(j,c+1):= Z(j,c); -- the real parts of the above 2 vecs are equal.
           end loop;
           if not (W_i(c+1) < Zero) then -- something went very wrong.
              raise Constraint_Error with "Fatal error in Unpack_Eig_Vectors";
           end if;
        end if;
     end loop;

  end Unpack_Eig_Vectors;

  ----------
  -- Norm --
  ----------

  function Norm
    (V_r, V_i     : in Col_Vector;
     Starting_Col : in Index := Index'First;
     Final_Col    : in Index := Index'Last)
     return Real
  is
     Max : Real := Zero;
     Sum, tst : Real := Zero;
     Scale  : Real := One;
     Scaled : Boolean := False;
     Scale_Exp : Integer;
     X, Y : Real;
  begin
     for i in Starting_Col .. Final_Col loop
        tst := Abs (V_r(i));
        if tst > Max then Max := tst; end if;
        tst := Abs (V_i(i));
        if tst > Max then Max := tst; end if;
     end loop;

     if Max > Sqrt_Max_Allowed_Real then
        Scale_Exp    := Integer'Min (Real'Exponent (Max), Real'Machine_Emax-4);
        Scale        := Two**(-Scale_Exp);
        Scaled       := True;
     end if;

     if Max < Sqrt_Min_Allowed_Real then
        Scale_Exp    := Integer'Max (Real'Exponent (Max), Real'Machine_Emin+4);
        Scale        := Two**(-Scale_Exp);
        Scaled       := True;
     end if;

     if Scaled then
        for k in Starting_Col .. Final_Col loop
           X   := Scale * V_r(k);
           Y   := Scale * V_i(k);
           Sum := Sum + X**2 + Y**2;
        end loop;
     else
        for k in Starting_Col .. Final_Col loop
           Sum := Sum + V_r(k)**2 + V_i(k)**2;
        end loop;
     end if;

     if Scaled then
        return Sqrt (Abs Sum) * Two**Scale_Exp;
     else
        return Sqrt (Abs Sum);
     end if;

  end Norm;

  ------------------------------
  -- Normalize_Column_Vectors --
  ------------------------------

  procedure Normalize_Column_Vectors
    (Z_r, Z_i : in out Matrix;
     Starting_Col : in Index;
     Final_Col    : in Index)
  is
     Min_Exp : constant Integer := Real'Machine_Emin/2 + Real'Machine_Emin/4;
     Min_Allowed_Real : constant Real := Two**Min_Exp;
     V_r, V_i : Col_Vector;
     Norm_Factor : Real;
  begin

     for Col_id in Starting_Col .. Final_Col loop

        for j in Starting_Col .. Final_Col  loop
           V_r(j) := Z_r(j, Col_id);
        end loop;

        for j in Starting_Col .. Final_Col  loop
           V_i(j) := Z_i(j, Col_id);
        end loop;

        Norm_Factor := 
           One / (Norm (V_r, V_i, Starting_Col, Final_Col) + Min_Allowed_Real);

        for j in Starting_Col .. Final_Col  loop
           Z_r(j, Col_id) := Z_r(j, Col_id) * Norm_Factor;
        end loop;

        for j in Starting_Col .. Final_Col  loop
           Z_i(j, Col_id) := Z_i(j, Col_id) * Norm_Factor;
        end loop;

     end loop;

  end Normalize_Column_Vectors;

  ----------------------------
  -- Apply_QR_to_Hessenberg --
  ----------------------------

  --  Apply_QR_to_Hessenberg is a translation of the algol procedure
  --  hqr2, num. math. 16, 181-204(1970) by Peters and Wilkinson.
  --  handbook for auto. comp., vol.ii-linear algebra,  372-395(1971).
  --
  --  Have  Z' A Z = H where H and Z are input below. (H is upper Hessenberg.) 
  --  The problem is to further elaborate Z so that  A Z = w Z.  The
  --  Col vectors of Z are eigvecs, and the diagonal elements of w are the
  --  eigvals.

  procedure Apply_QR_to_Hessenberg
    (H                    : in out Matrix;
     Z                    : in out Matrix;
     W_r, W_i             :    out Col_Vector;
     Norm                 :    out Real;
     Unscale_Eigs         :    out Real;
     Starting_Col         : in     Index   := Index'First;
     Final_Col            : in     Index   := Index'Last;
     Eigenvectors_Desired : in     Boolean := True;
     Id_of_Failed_Eig     :    out Integer)
  is
     N  : constant Real := Real(Final_Col) - Real(Starting_Col) + One;

     Max_Allowed_Iterations : constant Integer := 128 * Integer (N);

     Remaining_Iterations   : Integer := Max_Allowed_Iterations;
     Iteration_id           : Integer := 1;

     Emin             : constant Integer := Real'Machine_Emin;
     Min_Exp          : constant Integer := Emin - Emin / 16;
     Min_Allowed_Real : constant Real    := Two**(Min_Exp / 2 + 64);

     pragma Assert (Index'Base'First < Index'First - 1);
     pragma Assert (Final_Col > Starting_Col);

     Scale_Eigs : Real;
  begin

     W_r := (others => Zero); -- important init.
     W_i := (others => Zero);

     Id_of_Failed_Eig := Integer(Starting_Col) - 1; -- means none failed.

     --  compute matrix Norm

     Norm := Zero;
     for c in Starting_Col .. Final_Col loop
        for r in Starting_Col .. Final_Col loop
           Norm := Norm + Abs (H(r,c));
        end loop;
     end loop;

     Unscale_Eigs := One;
     if Norm > Two ** (Real'Machine_Emax / 2 - 2) then
        Scale_Eigs   := Two ** (- Real'Exponent (Norm));
        UnScale_Eigs := One / Scale_Eigs;
        for c in Starting_Col .. Final_Col loop
        for r in Starting_Col .. Final_Col loop
           H(r,c) := H(r,c) * Scale_Eigs;
        end loop;
        end loop;
     end if;

     Find_All_Eigenvalues:
     declare
        Eig_id_0  : Index'Base := Final_Col;
        Eig_id_2, Eig_id_1 : Index'Base;
        l_memo, m_memo : Index;
        Final_k : Boolean;
        t : Real := Zero; -- essential init
        tst1, tst2 : Real := Zero;
        Hmm : Real;
        x, y : Real;
        p,q,r,w,zz,xx : Real;
        Exp_s : Integer;
        Scale, Scale_Inv : Real;
        r_div_s, q_div_s, p_div_s : Real;
        q_div_ps, r_div_ps : Real;
        s, ss, tmp0 : Real := Zero;
     begin 

        Get_Next_Eigval:
        while Eig_id_0 >= Starting_Col loop

        -- iterate for eigenvalues:

        Iteration_id := 1;
        Eig_id_1     := Eig_id_0 - 1;
        Eig_id_2     := Eig_id_0 - 2;

        --  look for single small sub-diagonal element loop

        Find_Next_Root:
        loop

        for l in reverse Starting_Col .. Eig_id_0 loop

           l_memo := l;
           exit when l = Starting_Col;

           s := Abs (H(l-1,l-1)) + Abs (H(l,l));
           if s < Min_Allowed_Real then s := Norm; end if;

           exit when Arg_1_is_Negligible_Respect_Arg_2 (H(l,l-1), s);

        end loop;

        --  form shift.

        x := H(Eig_id_0,Eig_id_0);

        if l_memo = Eig_id_0 then                        -- 1 root found
           H(Eig_id_0,Eig_id_0) := x + t;
           W_r(Eig_id_0)        := H(Eig_id_0,Eig_id_0);
           W_i(Eig_id_0)        := Zero;
           Eig_id_0             := Eig_id_1;
           goto End_of_Loop_Get_Next_Eigval; 
           -- Notice we just decremented Eig_id_0.
        end if;

        y := H(Eig_id_1,Eig_id_1);
        w := H(Eig_id_0,Eig_id_1) * H(Eig_id_1,Eig_id_0);

        exit Find_Next_Root when l_memo = Eig_id_1;      -- 2 roots found

        if Remaining_Iterations = 0 then  -- didn't converge
           Id_of_Failed_Eig := Integer (Eig_id_0);
           return;
        end if;

        --  form exceptional shift

        if Iteration_id mod 14 = 0 then -- Forsythe_0 needs this.
           t := t + x;

           for i in Starting_Col .. Eig_id_0 loop
              H(i,i) := H(i,i) - x;
           end loop;

           s := Abs (H(Eig_id_0,Eig_id_1)) + Abs (H(Eig_id_1,Eig_id_2));
           x := 0.75 * s;
           y := x;
           w := -0.4375 * s * s;
        end if;

        Iteration_id := Iteration_id + 1;
        Remaining_Iterations := Remaining_Iterations - 1;

        --  look for two consecutive small loop sub-diagonal elements.

        for m in reverse l_memo .. Eig_id_2 loop
           m_memo := m;
           Hmm := H(m,m);
           r := x - Hmm;
           s := y - Hmm;
           p := (r * s - w) / H(m+1,m) + H(m,m+1);
           q := H(m+1,m+1) - Hmm - r - s;
           r := H(m+2,m+1);

           Scale := One / (Abs (p) + Abs (q) + Abs (r) + Min_Allowed_Real);
           p := p * Scale;
           q := q * Scale;
           r := r * Scale;

           exit when m = l_memo;

           tst1 := Abs (H(m,m-1)) * (Abs (q) + Abs (r));
           tst2 := Abs (p)*(Abs (H(m-1,m-1)) + Abs (Hmm) + Abs (H(m+1,m+1)));

           exit when Arg_1_is_Negligible_Respect_Arg_2 (tst1, tst2);

        end loop;

        for i in m_memo+2 .. Eig_id_0 loop
           H(i,i-2) := Zero;
           if (i /= m_memo+2) then
              H(i,i-3) := Zero;
           end if;
        end loop;

        -- double qr step w/ rows  l_memo .. Eig_id_0,  columns  m_memo .. Eig_id_0

        Double_QR:
        for k in m_memo .. Eig_id_1 loop
           Final_k := (k = Eig_id_1);

           if k /= m_memo then

              p := H(k,k-1);
              q := H(k+1,k-1);
              if Final_k then
                 r := Zero;
              else
                 r := H(k+2,k-1);
              end if;

              Scale := Abs (p) + Abs (q) + Abs (r);
              if Scale < Min_Allowed_Real then
                 goto End_of_Loop_Double_QR;   -- and then increment k.
              end if;

              if Scale < Two**(Emin / 2)  or else  Scale > Two**(Real'Emax / 4) then
                 Exp_s     := Real'Exponent (Scale);
                 scale     := Two ** (Exp_s);
                 Scale_Inv := Two ** (-Exp_s);
                 p := p * Scale_Inv;
                 q := q * Scale_Inv;
                 r := r * Scale_Inv;
              else 
                 scale     := One;
                 scale_inv := One;
              end if;

           end if;

           ss := p*p + q*q + r*r;
           s  := Real'Copy_Sign (Sqrt (ss), p);

           if (k /= m_memo) then
              H(k,k-1) := -s * Scale;
           else
              if (l_memo /= m_memo) then
                 H(k,k-1) := -H(k,k-1);
              end if;
           end if;

           p_div_s  := p / s;
           q_div_s  := q / s;
           r_div_s  := r / s;
           q_div_ps := q / (p + s);
           r_div_ps := r / (p + s);

           if not Final_k then

              --  rows

              for j in k .. Final_Col loop
                 tmp0     :=  H(k,j) + (q_div_ps * H(k+1,j) + r_div_ps * H(k+2,j));
                 H(k,j)   := -p_div_s * H(k,j) - (q_div_s*H(k+1,j) + r_div_s*H(k+2,j));
                 H(k+1,j) :=  H(k+1,j) - tmp0 * q_div_s;
                 H(k+2,j) :=  H(k+2,j) - tmp0 * r_div_s;
              end loop;

              --  columns

              for i in Starting_Col .. Index'Min (Eig_id_0,k+3) loop
                 tmp0     :=  H(i,k) + (q_div_ps * H(i,k+1) + r_div_ps * H(i,k+2));
                 H(i,k)   := -p_div_s * H(i,k) - (q_div_s*H(i,k+1) + r_div_s*H(i,k+2));
                 H(i,k+1) :=  H(i,k+1) - tmp0 * q_div_s;
                 H(i,k+2) :=  H(i,k+2) - tmp0 * r_div_s;
              end loop;

              --  accumulate transformations

              if Eigenvectors_Desired then
              for i in Starting_Col .. Final_Col loop
                 tmp0     :=  Z(i,k) + (q_div_ps * Z(i,k+1) + r_div_ps * Z(i,k+2));
                 Z(i,k)   := -p_div_s * Z(i,k) - (q_div_s*Z(i,k+1) + r_div_s*Z(i,k+2));
                 Z(i,k+1) :=  Z(i,k+1) - tmp0 * q_div_s;
                 Z(i,k+2) :=  Z(i,k+2) - tmp0 * r_div_s;
              end loop;
              end if;

           else

              for j in k .. Final_Col loop
                 tmp0     :=  H(k,j) + q_div_ps * H(k+1,j);
                 H(k,j)   := -p_div_s * H(k,j) - q_div_s * H(k+1,j);
                 H(k+1,j) :=  H(k+1,j) - q_div_s * tmp0;
              end loop;

              --  columns

              for i in Starting_Col .. Index'Min (Eig_id_0,k+3) loop
                 tmp0     :=  H(i,k) + q_div_ps * H(i,k+1);
                 H(i,k)   := -p_div_s * H(i,k) - q_div_s * H(i,k+1);
                 H(i,k+1) :=  H(i,k+1) - q_div_s * tmp0;
              end loop;

              --  accumulate transformations

              if Eigenvectors_Desired then
              for i in Starting_Col .. Final_Col loop
                 tmp0     :=  Z(i,k) + q_div_ps * Z(i,k+1);
                 Z(i,k)   := -p_div_s * Z(i,k) - q_div_s * Z(i,k+1);
                 Z(i,k+1) :=  Z(i,k+1) - q_div_s * tmp0;
              end loop;
              end if;

           end if;

           <<End_of_Loop_Double_QR>>  null;

        end loop Double_QR; -- increment k

        end loop Find_Next_Root;

        --  two roots found. x, y, w were defined up top as:
        --
        -- x := H(Eig_id_0,Eig_id_0);
        -- y := H(Eig_id_1,Eig_id_1);
        -- w := H(Eig_id_0,Eig_id_1) * H(Eig_id_1,Eig_id_0);

        p  := Half * (y - x);
        q  := p * p + w;
        zz := Sqrt (Abs (q));

        H(Eig_id_0,Eig_id_0) := x + t;
        H(Eig_id_1,Eig_id_1) := y + t;
        x                    := x + t;

        if q >= Zero then          --  real pair
         
           zz            := p + Real'Copy_Sign (zz, p);
           W_r(Eig_id_1) := x + zz;
           W_r(Eig_id_0) := W_r(Eig_id_1);
           if (Abs zz > Min_Allowed_Real) then
              W_r(Eig_id_0) := x - w / zz;
           end if;
           W_i(Eig_id_1) := Zero;
           W_i(Eig_id_0) := Zero;

           xx := H(Eig_id_0,Eig_id_1);

           declare 
              hypot, sn_lo, sn_hi, cs_lo, cs_hi : Real;
              tmp1, tmp0 : Real;
           begin
         
              Rotation_Factors_2 (xx, zz, sn_lo, sn_hi, cs_lo, cs_hi, hypot);
              -- sn = a / Sqrt(a*a + b*b), cs = b / Sqrt(a*a + b*b), (a=xx, b = zz)
         
              -- rows
              for j in Eig_id_1 .. Final_Col loop
                 tmp1          := H(Eig_id_1,j);
                 tmp0          := H(Eig_id_0,j);
                 H(Eig_id_1,j) := sn_hi*tmp0 + cs_hi*tmp1 + (cs_lo*tmp1 + sn_lo*tmp0);
                 H(Eig_id_0,j) :=-sn_hi*tmp1 + cs_hi*tmp0 + (cs_lo*tmp0 - sn_lo*tmp1);
              end loop;
         
              --  columns
              for i in Starting_Col .. Eig_id_0 loop
                 tmp1          := H(i,Eig_id_1);
                 tmp0          := H(i,Eig_id_0);
                 H(i,Eig_id_1) := sn_hi*tmp0 + cs_hi*tmp1 + (cs_lo*tmp1 + sn_lo*tmp0);
                 H(i,Eig_id_0) :=-sn_hi*tmp1 + cs_hi*tmp0 + (cs_lo*tmp0 - sn_lo*tmp1);
              end loop;
         
              --  columns
              for i in Starting_Col .. Final_Col loop
                 tmp1          := Z(i,Eig_id_1);
                 tmp0          := Z(i,Eig_id_0);
                 Z(i,Eig_id_1) := sn_hi*tmp0 + cs_hi*tmp1 + (cs_lo*tmp1 + sn_lo*tmp0);
                 Z(i,Eig_id_0) :=-sn_hi*tmp1 + cs_hi*tmp0 + (cs_lo*tmp0 - sn_lo*tmp1);
              end loop;
           end;

        else   --  complex pair
          
            W_r(Eig_id_1) := p + x;
            W_r(Eig_id_0) := p + x;
            W_i(Eig_id_1) := zz;
            W_i(Eig_id_0) :=-zz;

        end if;

        Eig_id_0 := Eig_id_2;

        <<End_of_Loop_Get_Next_Eigval>> null;

        end loop Get_Next_Eigval;

     end Find_All_Eigenvalues;

  end Apply_QR_to_Hessenberg;

  ---------------------------
  -- Find_all_Eigenvectors --
  ---------------------------

  --  Backsubstitute to find vectors of upper triangular form H.

  procedure Find_all_Eigenvectors
    (H                    : in out Matrix;
     Z                    : in out Matrix;
     Z_i                  :    out Matrix;
     W_r, W_i             : in     Col_Vector;
     Norm                 : in     Real;
     Starting_Col         : in     Index   := Index'First;
     Final_Col            : in     Index   := Index'Last)
  is
     Eig_Start, Eig_Memo : Index;
     Eig_2, Eig_1 : Index'Base;
     ra, sa : Real;
     p, q, r, s, w, zz : Real;
     Sum : Real;

     Min_Exp          : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real    := Two**(Min_Exp - Min_Exp/16);

     pragma Assert (Index'Base'First < Index'First - 1);
     pragma Assert (Final_Col > Starting_Col);
  begin
     if Norm = Zero then
        Z_i := (others => (others => Zero));
        return;
     end if;

     Get_Next_Eigenvector:
     for Eig_0 in reverse Starting_Col .. Final_Col loop

        p := W_r(Eig_0);
        q := W_i(Eig_0);

        if Abs W_i(Eig_0) = Zero then  -- Make real vector

           H(Eig_0,Eig_0) := One;
           Eig_Memo       := Eig_0;

           if Eig_0 > Starting_Col then

           Get_Real_Vector:
           for i in reverse Starting_Col .. Eig_0-1 loop
              w := H(i,i) - p;

              r := Zero;
              for j in Eig_Memo .. Eig_0 loop
                 r := r + H(i,j) * H(j,Eig_0);
              end loop;

              if W_i(i) < Zero then 

                 zz := w;
                 s  := r;

              else

                 Eig_Memo := i;
                 if (Abs W_i(i) < Min_Allowed_Real) then
                    declare
                       t : Real := w;
                       Exp_t : Integer;
                    begin
                       if Abs t < Min_Allowed_Real then
                          Exp_t := Real'Machine_Mantissa + 27;
                          t  := Norm * Two**(-Exp_t) + Min_Allowed_Real;
                       end if;
                       H(i,Eig_0) := -r / t;
                    end;

                 else

                    --  solve real equations

                    declare
                       x, y, q, t : Real;
                    begin
                       x := H(i,i+1);
                       y := H(i+1,i);
                       q := (W_r(i) - p)**2 + W_i(i)**2;
                       t := (x * s - zz * r) / q;
                       H(i,Eig_0) := t;
                       if (Abs (x) > Abs (zz)) then
                          H(i+1,Eig_0) := (-r - w * t) / x;
                       else
                          H(i+1,Eig_0) := (-s - y * t) / zz;
                       end if;
                    end;

                 end if;

                 --  overflow control
                 declare
                    t, Scale : Real;
                    Exp_t : Integer;
                 begin
                    t := Abs (H(i,Eig_0));
                    if t > Min_Allowed_Real then  -- "t > Zero" fails
                       Exp_t := Real'Exponent (t);
                       if Exp_t > Real'Machine_Mantissa + 27 then
                          Scale := Two**(-Exp_t);
                          for j in i .. Eig_0 loop
                             H(j,Eig_0) := H(j,Eig_0) * Scale;
                          end loop;
                       end if;
                    end if;
                 end;

              end if;
           end loop Get_Real_Vector;  -- in i
           end if;                    -- if Eig_0 > Starting_Col

           -- Now go to end of loop Get_Next_Eigenvector, and increment Eig_0.

        elsif q < Zero then

           --  Complex vectors. Last vector component chosen imaginary
           --  so that eigenvector matrix is triangular.

           Eig_1 := Eig_0 - 1;
           Eig_2 := Eig_0 - 2;

           -- check

           if (Abs (H(Eig_0,Eig_1)) > Abs (H(Eig_1,Eig_0))) then
              H(Eig_1,Eig_1) := q / H(Eig_0,Eig_1);
              H(Eig_1,Eig_0) := -(H(Eig_0,Eig_0) - p) / H(Eig_0,Eig_1);
           else
              Div (Zero, -H(Eig_1,Eig_0),
                   H(Eig_1,Eig_1)-p, q,
                   H(Eig_1,Eig_1), H(Eig_1,Eig_0));
           end if;
           H(Eig_0,Eig_1) := Zero;
           H(Eig_0,Eig_0) := One;
           Eig_Start      := Eig_1;

           Get_Complex_Vector:
           for i in reverse Starting_Col .. Eig_2 loop

              w  := H(i,i) - p;

              ra := Zero;
              sa := Zero;
              for j in Eig_Start .. Eig_0 loop
                 ra := ra + H(i,j) * H(j,Eig_1);
                 sa := sa + H(i,j) * H(j,Eig_0);
              end loop;

              if (W_i(i) < Zero) then
                 zz := w;
                 r  := ra;
                 s  := sa;
                 goto End_of_Loop_Get_Complex_Vector; 
              end if;

              Eig_Start := i; -- decrement start of previous loop

              if (Abs W_i(i) < Min_Allowed_Real) then

                 Div (-ra, -sa, w, q, H(i,Eig_1), H(i,Eig_0));

              else

                 --  solve complex equations

                 declare
                    x, y, vr, vi, tst1 : Real;
                 begin
                    x := H(i,i+1);
                    y := H(i+1,i);
                    vr := (W_r(i) - p) * (W_r(i) - p) + W_i(i) * W_i(i) - q * q;
                    vi := (W_r(i) - p) * Two * q;
                    if (Abs vr < Min_Allowed_Real and
                        Abs vi < Min_Allowed_Real)
                    then
                       tst1 := Norm * (Abs(w)+Abs(q)+Abs(x)+Abs(y)+Abs(zz));
                       vr   := tst1 * (Real'Epsilon * Two**(-10));
                    end if;

                    Div(x*r-zz*ra+q*sa, x*s-zz*sa-q*ra,
                        vr, vi,
                        H(i,Eig_1), H(i,Eig_0));

                    if (Abs (x) > Abs (zz) + Abs (q)) then
                       H(i+1,Eig_1) := (-ra - w * H(i,Eig_1) + q * H(i,Eig_0)) / x;
                       H(i+1,Eig_0) := (-sa - w * H(i,Eig_0) - q * H(i,Eig_1)) / x;
                    else
                       Div(-r-y*H(i,Eig_1), -s-y*H(i,Eig_0),
                           zz, q,
                           H(i+1,Eig_1), H(i+1,Eig_0));
                    end if;
                 end;
              end if;

              --  overflow control
              declare
                 t, Scale : Real;
                 Exp_t : Integer;
              begin
                 t := Real'Max (Abs (H(i,Eig_1)), Abs (H(i,Eig_0)));
                 if t > Min_Allowed_Real then
                    Exp_t := Real'Exponent (t);
                    if Exp_t > Real'Machine_Mantissa + 27 then
                       Scale := Two**(-Exp_t);
                       for j in i .. Eig_0 loop
                          H(j,Eig_1) := H(j,Eig_1) * Scale;
                          H(j,Eig_0) := H(j,Eig_0) * Scale;
                       end loop;
                    end if;
                 end if;
              end;

              <<End_of_Loop_Get_Complex_Vector>> null;

        end loop Get_Complex_Vector;       -- increment i

     end if; -- Abs W_i(Eig_0) < Min_Allowed_Real => do real vecs, else complex.

     end loop Get_Next_Eigenvector;        -- increment Eig_0

     --  multiply by transformation matrix to get the eigenvectors.
     --  Have  H v = w v, where Z' A Z = H.  Then Z H Z'(Zv) = w (Zv)
     --  implies that Z*v are the eigvecs of A = Z H Z'.  (To save much
     --  space, the "v" eigenvectors were stored in H). 
   
     --for j in reverse Starting_Col .. Final_Col loop
     --   for i in Starting_Col .. Final_Col loop
     --      Sum := Zero;
     --      for k in Starting_Col .. j loop
     --         Sum := Sum + Z(i,k) * H(k,j);
     --      end loop;
     --      Z(i,j) := Sum;
     --   end loop;
     --end loop;
 
     declare 
        Z_row_vec, Z_new_row : array (Starting_Col .. Final_Col) of Real; 
     begin
        for r in Starting_Col .. Final_Col loop
           for k in Starting_Col .. Final_Col loop  -- reuse a row of Z below
              Z_row_vec(k) :=  Z(r, k);
           end loop;
           for c in Starting_Col .. Final_Col loop
              Sum := Zero;
              for k in Starting_Col .. c loop       -- H is upper triangular
               --Sum := Sum + Z(r, k) * H(k, c);
                 Sum := Sum + Z_row_vec(k) * H(k, c); -- H is Convention(Fortran).
              end loop;
              Z_new_row(c) := Sum;
           end loop;
           for c in Starting_Col .. Final_Col loop
              Z(r, c) :=  Z_new_row(c);
           end loop;
        end loop;
     end;

     Unpack_Eig_Vectors (W_i, Z, Z_i, Starting_Col, Final_Col);
     -- Gets (Z_r, Z_i) from Z. Destroys old Z.
     -- The out parameter Z_i is finally initialized here.

  end Find_all_Eigenvectors;

  ------------------------
  -- Sort_Eigs_And_Vecs --
  ------------------------

  procedure Sort_Eigs_And_Vecs
    (Z_r, Z_i     : in out Matrix;       -- eigvecs are columns of Z
     W_r, W_i     : in out Col_Vector;   -- eigvals are vectors W
     Starting_Col : in     Index   := Index'First;
     Final_Col    : in     Index   := Index'Last)
  is
     Max_Eig, tmp : Real;
     Eig_Size : Col_Vector := (others => Zero);
     i_Max : Index;
  begin
     if Starting_Col < Final_Col then

     for i in Starting_Col .. Final_Col loop
        Eig_Size(i) := Hypo.Hypotenuse (W_r(i), W_i(i));
     end loop;

     for i in Starting_Col .. Final_Col-1 loop

        Max_Eig := Eig_Size(i);
        i_Max   := i;

        for j in i+1 .. Final_Col loop
           if Eig_Size(j) > Max_Eig then
              Max_Eig := Eig_Size(j);
              i_Max   := j;
           end if;
        end loop;

        if i_Max > i then

           tmp        := W_r(i);
           W_r(i)     := W_r(i_Max);
           W_r(i_Max) := tmp;

           tmp        := W_i(i);
           W_i(i)     := W_i(i_Max);
           W_i(i_Max) := tmp;

           tmp             := Eig_Size(i);
           Eig_Size(i)     := Eig_Size(i_Max);
           Eig_Size(i_Max) := tmp;

           -- swap cols of Z, the eigenvectors:

           for k in Starting_Col .. Final_Col loop
              tmp           := Z_r(k, i);
              Z_r(k, i)     := Z_r(k, i_Max);
              Z_r(k, i_Max) := tmp;
           end loop;

           for k in Starting_Col .. Final_Col loop
              tmp           := Z_i(k, i);
              Z_i(k, i)     := Z_i(k, i_Max);
              Z_i(k, i_Max) := tmp;
           end loop;

        end if;

     end loop;
     end if;

  end Sort_Eigs_And_Vecs;

  procedure Transpose 
    (M            : in out Matrix;
     Starting_Col : in     Index        := Index'First;
     Final_Col    : in     Index        := Index'Last)
  is
     tmp : Real;
  begin
     for r in Starting_Col .. Final_Col loop
        for c in Starting_Col .. r loop
           tmp    := M(c,r);
           M(c,r) := M(r,c);
           M(r,c) := tmp;
        end loop;
     end loop;
  end Transpose;

  ---------------
  -- Decompose --
  ---------------

  procedure Decompose
    (A                    : in out Matrix;
     Z_r, Z_i             :    out Matrix;
     W_r, W_i             :    out Col_Vector;
     Id_of_Failed_Eig     :    out Integer;
     Starting_Col         : in     Index        := Index'First;
     Final_Col            : in     Index        := Index'Last;
     Eigenvectors_Desired : in     Boolean      := True;
     Balance_Policy       : in     Balance_Code := Disabled)
  is
     Norm : Real := One;
     Unscale_Eigs : Real := One;
     D : Col_Vector := (others => One);
  begin
     Id_of_Failed_Eig := Integer (Final_Col); -- init out parameter.

     Z_r := Givens_Hess.Identity;

     if Final_Col <= Starting_Col then
        raise Constraint_Error with "Can't have Final_Col <= Starting_Col";
     end if;

     Balance_Matrix (A, D, Starting_Col, Final_Col, Balance_Policy);

     -- if matrix was balanced, then multiply D by I to get starting
     -- point for calculation of the Q matrix: Initial_Q => Z_r.
     -- 
     -- ** The Q matrices won't be orthogonal if the matrix is balanced. **

     for c in Index loop
        Z_r(c,c) := D(c);
     end loop;

     -- More generally, if  Initial_Q is not I, but some Z_r:
     -- 
     --for r in Starting_Col .. Final_Col loop
     --for c in Starting_Col .. Final_Col loop
        --Z_r(r,c) := Z_r(r,c) * D(c);
     --end loop;
     --end loop;

     Givens_Hess.Upper_Hessenberg
       (A            => A,
        Q            => Z_r,
        Starting_Col => Starting_Col,
        Final_Col    => Final_Col,
        Initial_Q    => Z_r);

     -- A is now in Upper_Hessenberg form, thanks to Z_r: A_true = Z_r' A Z_r.
     -- Procedure Upper_Hessenberg initializes "out" parameter Z_r.

     for i in 1 .. 8 loop
        Initial_QR.Lower_Diagonal_QR_Iteration
          (A            => A,
           Q            => Z_r,
           Shift        => Zero,
           Starting_Col => Starting_Col,
           Final_Col    => Final_Col);
     end loop;
     -- Helps eigvec calculation with several matrices, (eg. Pas_Fib, Vandermonde)

     -- Next, Z_r is input to Apply_QR_to_Hessenberg, where it is (optionally)
     -- transformed into the set of eigenvectors of A.
     --
     -- A is now in Upper_Hessenberg form, thanks to Z_r: A_true = Z_r A Z_r'.

     Apply_QR_to_Hessenberg
       (A,
        Z_r, 
        W_r, W_i,
        Norm,
        Unscale_Eigs,
        Starting_Col, Final_Col,
        Eigenvectors_Desired,     -- if true then Z_i will be initialized
        Id_of_Failed_Eig);

     -- If  Eigenvectors_Desired = False  then  Z_i remains uninitialized,
     -- (which may free up space with large matrices).

     if Eigenvectors_Desired then

        Find_all_Eigenvectors
          (A,
           Z_r, Z_i,
           W_r, W_i,
           Norm,
           Starting_Col, Final_Col);

        Normalize_Column_Vectors (Z_r, Z_i, Starting_Col, Final_Col);

        for c in Starting_Col .. Final_Col loop
           W_r(c) := Unscale_Eigs * W_r(c);
           W_i(c) := Unscale_Eigs * W_i(c);
        end loop;

     end if;

  end Decompose;

end Peters_Eigen;
