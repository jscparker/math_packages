
-----------------------------------------------------------------------
-- package body e_Jacobi_Eigen, extended precision Jacobi eigen-decomposition
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


package body e_Jacobi_Eigen is

  Reciprocal_Epsilon : constant Real := One / e_Real_Model_Epsilon;

  ---------------------------------
  -- Get_Jacobi_Rotation_Factors --
  ---------------------------------

  --  all underflows are OK here.
  --  no   overflows are OK here.
  --  so we test for  Q / P  overflows by calculating  P / Q

  procedure Get_Jacobi_Rotation_Factors
    (P, Q    : in     Real;
     s       :    out Real;
     tau     :    out Real;
     Delta_D :    out Real)
  is 
     t, Gamma : Real;  
     -- function "/" (x, y : Real) return Real renames Divide;
     -- faster than stnd "/", but scarcely matters.
  begin

     s       := Zero;
     tau     := Zero;
     Delta_D := Zero;

     if Abs (Q) > e_Real_Safe_Min and then
        Abs (P) > e_Real_Safe_Min and then
        Abs (Q) < e_Real_Safe_Max and then
        Abs (P) < e_Real_Safe_Max then

        --  Following Gamma is usually the more accurate.

        Gamma := P / (Abs (Q) + Sqrt (P*P + Q*Q));

	if Q < Zero then Gamma := -Gamma; end if;

        -- usual case overwhelmingly. If you scale matrix to unit Norm,
	-- then no overflows because Matrix norm is preserved, preventing
	-- large P, Q if they are not large to begin with. (So the above
	-- tests for < e_Real_Safe_Max would be unnecessary.)
	-- (Requires scaling of matrix prior to decomposition to be
	-- absolutely sure.)

        -- Should be able to do the following w/o any tests of p,q if don't care
	-- about quality of answer in P < e_Safe_Small, Q < e_Safe_Small limit. 
	--
        --Gamma := P / (Abs(Q) + Sqrt (P*P + Q*Q) + e_Safe_Small);
	--if Q < Zero then Gamma := -Gamma; end if;

     elsif  Abs (Q) > Abs (P) then 
     
        -- Abs (P) > 0 was a tested before arrival, which implies Abs (Q) > 0.

        t       := P / Q;
        Gamma   := t / (One + Sqrt (One + t*t));
        --  Underflow OK; overflow not allowed.

     elsif Abs (P) >= Abs (Q) then 

        t       := Q / P;
        Gamma   := One / (Abs (t) + Sqrt (One + t*t));
        --  Underflow OK; overflow not allowed.

	if t < Zero then Gamma := -Gamma; end if;

     else

        return;  -- must have hit some inf's. Use stnd rotation init'd above.

     end if;

     declare
        c     : Real;
     begin
        c       := Reciprocal_Sqrt (One + Gamma*Gamma); -- Cosine (ok)
      --c       := Sqrt (One / (One + Gamma*Gamma));    -- Cosine
        s       := c * Gamma;                           -- Sine
        tau     := s / (One + c);                       -- -cos_minus_1_over_sin
        Delta_D := P * Gamma;
     end;

  end Get_Jacobi_Rotation_Factors;

  ---------------------
  -- Eigen_Decompose --
  ---------------------

  procedure Eigen_Decompose
    (A : in out Matrix;
     Q_tr : out Matrix;
     Eigenvals              : out Col_Vector;
     No_of_Sweeps_Performed : out Natural;
     Total_No_of_Rotations  : out Natural;
     Start_Col              : in Index := Index'First;
     Final_Col              : in Index := Index'Last;
     Eigenvectors_Desired   : in Boolean      := False)
  is
     D  : Col_Vector renames Eigenvals;
     Z, B  : Col_Vector;

     Max_Allowed_No_of_Sweeps : constant Positive := 256; -- badly scaled need lots
     No_of_Preliminary_Sweeps : constant Positive := 14;

     --Reciprocal_Epsilon : constant Real := One / e_Real_Model_Epsilon;
     -- Good stnd setting for the effective eps is:  Real'Epsilon * Two**(-2).
     -- Usually, Real'Epsilon := 2.0**(-50) for 15 digit Reals.

     Matrix_Size : constant Real_8 := Real_8(Final_Col) - Real_8(Start_Col) + 1.0;
     No_of_Off_Diag_Elements : constant Real_8 := 0.5*Matrix_Size*(Matrix_Size-1.0);
     Mean_Off_Diagonal_Element_Size : Real;

     s, g, h, tau : Real;   -- Rutishauser variable names.
     Q, Delta_D   : Real;
     Sum, Pivot, Threshold : Real;
  begin

     -- Initialize all out parameters. D renames Eigenvals.
     -- Q_tr starts as Identity; is rotated into set of Eigenvectors of A.
 
     Q_tr := (others => (others => Zero));
     for j in Index loop
        Q_tr(j,j) := One;
     end loop;
 
     Z := (others => Zero);
     B := (others => Zero);
     D := (others => Zero);
     for j in Start_Col .. Final_Col loop  -- assume A not all init
        D(j) := A(j,j);
        B(j) := A(j,j);
     end loop;
 
     No_of_Sweeps_Performed := 0;
     Total_No_of_Rotations  := 0;
 
     if Matrix_Size <= 1.0 then  return;   end if; -- right answer for Size=1.
 
 
     Sweep: for Sweep_id in 1 .. Max_Allowed_No_of_Sweeps loop
 
       Sum := Zero;
       for N in Start_Col .. Final_Col-1 loop    --sum off-diagonal elements
       for I in N+1 .. Final_Col loop
          Sum := Sum + Abs (A(N,I));
       end loop;
       end loop;
       Mean_Off_Diagonal_Element_Size := Sum / (+No_of_Off_Diag_Elements);


       exit Sweep when Mean_Off_Diagonal_Element_Size < Min_Allowed_Real;

 
       if Sweep_id > No_of_Preliminary_Sweeps then
          Threshold := Zero;
       else
          Threshold := One * Mean_Off_Diagonal_Element_Size;
       end if;
 
       for N in Start_Col .. Final_Col-1 loop
       for I in N+1 .. Final_Col loop
 
         Pivot := A(N,I);
 
         --  Have to zero out sufficiently small A(I,N) to get convergence,
	 --  ie, to get  Off_Diag_Sum -> 0.0.
	 --  After 4 sweeps all A(I,N) are small so that 
	 --  A(I,N) / Epsilon  will never overflow. The test is
	 --  A(I,N) / Epsilon <= Abs D(I) and A(I,N) / Epsilon <= Abs D(N).
 
         if
	    (Sweep_id >          No_of_Preliminary_Sweeps) and then
            (Reciprocal_Epsilon * Abs (Pivot) <= Abs D(N)) and then
            (Reciprocal_Epsilon * Abs (Pivot) <= Abs D(I))
         then
 
            A(N,I) := Zero;
 
         elsif Abs (Pivot) > Threshold then
 
            Q := Half * (D(I) - D(N));
  
            Get_Jacobi_Rotation_Factors (Pivot, Q, s, tau, Delta_D);
            -- Pivot=A(N,I)
  
            D(N) := D(N) - Delta_D;  -- Locally D is only used for threshold test.
            D(I) := D(I) + Delta_D;
            Z(N) := Z(N) - Delta_D;  -- Z is reinitialized to 0 each sweep, so Z
            Z(I) := Z(I) + Delta_D;  -- sums the small d's 1st. Helps a tad.
  
            A(N,I) := Zero;
  
            for j in Start_Col .. N-1 loop
               g := A(j,N);
               h := A(j,I);
               A(j,N) := g-s*(h+g*tau);
               A(j,I) := h+s*(g-h*tau);
            end loop;
            for j in N+1 .. I-1 loop
               g := A(N,j);
               h := A(j,I);
               A(N,j) := g-s*(h+g*tau);
               A(j,I) := h+s*(g-h*tau);
            end loop;
            for j in I+1 .. Final_Col loop
               g := A(N,j);
               h := A(I,j);
               A(N,j) := g-s*(h+g*tau);
               A(I,j) := h+s*(g-h*tau);
            end loop;
  
            if Eigenvectors_Desired then
            for j in Start_Col .. Final_Col loop
               g := Q_tr(N,j);
               h := Q_tr(I,j);
               Q_tr(N,j) := g-s*(h+g*tau);
               Q_tr(I,j) := h+s*(g-h*tau);
            end loop;
            end if;
  
            Total_No_of_Rotations := Total_No_of_Rotations + 1;
 
         end if; -- if (Sweep_id > No_of_Preliminary_Sweeps)
 
       end loop; --I loop (Col)
       end loop; --N loop (Row)
 
       for j in Start_Col .. Final_Col loop  -- assume A not all initialized
          B(j) := B(j) + Z(j);
          D(j) := B(j);
          Z(j) := Zero;
       end loop;
 
     end loop Sweep; --Sweep_id loop
 
  end Eigen_Decompose;

  ---------------
  -- Sort_Eigs --
  ---------------

  procedure Sort_Eigs
    (Eigenvals         : in out Col_Vector;
     Q_tr              : in out Matrix;             -- rows  are the eigvectors
     Start_Col         : in     Index   := Index'First;
     Final_Col         : in     Index   := Index'Last;
     Sort_Eigvecs_Also : in     Boolean := False)
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
           tmp            := Q_tr(i,k);
           Q_tr(i,k)      := Q_tr(Max_id,k);
           Q_tr(Max_id,k) := tmp;
        end loop;
      end if;

    end loop;
    end if;

  end Sort_Eigs;

end e_Jacobi_Eigen;

