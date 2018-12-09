
-----------------------------------------------------------------------
-- package body Spline, cubic splines
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

package body Spline is

  Zero      : constant Real := +0.0;
  One       : constant Real := +1.0;
  Two       : constant Real := +2.0;
  Three     : constant Real := +3.0;
  Six       : constant Real := +6.0;
  One_Third : constant Real := One / Three;

  --------------------
  -- Prepare_X_Data --
  --------------------
 
  procedure Prepare_X_Data 
    (X_Stuff     :    out X_Structure;
     X_Data      : in     Data_Vector;
     I_Start     : in     Index    := Index'first;
     I_Finish    : in     Index    := Index'Last;
     Bound_First : in     Boundary := Natural_BC;
     Bound_Last  : in     Boundary := Natural_BC)
  is
     D_1st  : constant Tri.DiagonalID := Tri.DiagonalID'First;
     D_0    : constant Tri.DiagonalID := Tri.DiagonalID'Succ(Tri.DiagonalID'First);
     D_Last : constant Tri.DiagonalID := Tri.DiagonalID'Last;
  
     No_Of_Points_Minus_1 : constant Index'Base := I_Finish - I_Start;
  
     S : X_Structure renames X_Stuff;
  begin
 
     if No_Of_Points_Minus_1 < 3 then
        raise Constraint_Error with "Too few points for spline calculation";
     end if;
   
     for i in I_Start .. I_Finish-1 loop
        S.dX(i) := X_Data(i+1) - X_Data(i);
        if S.dX(i) < Zero then
           raise Constraint_Error with "X axis points are out of order.";
        end if;
        if Abs (S.dX(i)) < Two ** (Real'Machine_Emin / 2) then
           raise Constraint_Error with "Points too close together on X axis.";
        end if;
     end loop;
  
     for i in I_Start .. I_Finish-1 loop
        S.dX_Inverse(i) := One / S.dX(i);
     end loop;
  
     -- We need to solve the matrix equation M*Z = V, for the vector Z.  
     -- First initialize M and V.  M is Tridiagonal. The lower diagonal
     -- is M(-1, .. ). The elements of the diagonal "vectors" are
     -- indexed by the row index the element would have in
     -- matrix form, so M(-1,*) starts at I_Start+1, etc.
     --
     -- Impose a value Boundary_Value_First = Alpha*Y_dot_dot + Beta*Y_dot:
     --
     --  (2*h_0*Beta - 6*Alpha)*C_0 + Beta*h_0*C_1
     --         = -3*Boundary_Value_First + 3*Beta*(Y_1 - Y_0)/h_0.
     --
     -- At the other end set Boundary_Value_Last = Alpha2*Y_dot_dot + Beta2*Y_dot,
     -- to get:
     --
     --  (2*h_n-1*Beta2 + 6*Alpha2)*C_n + Beta2*h_n-1*C_n-1
     --         = 3*Boundary_Value_Last - 3*Beta2*(Y_n - Y_n-1)/h_n-1.
    
     S.M(D_Last)(I_Start) := S.dX(I_Start) * Bound_First.Beta;
     S.M(D_0)(I_Start)    := S.dX(I_Start) * Bound_First.Beta * Two
            - Six * Bound_First.Alpha;
  
     S.M(D_1st)(I_Finish) := S.dX(I_Finish-1) * Bound_Last.Beta;
     S.M(D_0)(I_Finish)   := S.dX(I_Finish-1) * Bound_Last.Beta * Two
             + Six * Bound_Last.Alpha;
  
     -- Get rest of matrix M and vector V.
     -- Clamped and Natural_BC Spline types have the following in common.
  
     for i in I_Start+1 .. I_Finish-1 loop
        S.M(D_Last)(i) := S.dX(i);
     end loop;
  
     for i in I_Start+1 .. I_Finish-1 loop
        S.M(D_1st)(i) := S.dX(i-1);
     end loop;
  
     for i in I_Start+1 .. I_Finish-1 loop
        S.M(D_0)(i) := Two * (S.dX(i) + S.dX(i-1));
     end loop;
  
     Tri.LU_Decompose (S.M, I_Start, I_Finish);
 
     S.I_Start     := I_Start;
     S.I_Finish    := I_Finish;
     S.Bound_First := Bound_First;
     S.Bound_Last  := Bound_Last;
     S.Initialized := True;
 
  end Prepare_X_Data;
 
  ----------------
  -- Get_Spline --
  ----------------
 
  procedure Get_Spline 
    (Spline  :    out Spline_Coefficients;
     X_Stuff : in     X_Structure;
     Y_Data  : in     Data_Vector)
  is
     dY         : Data_Vector    := (others => Zero);
     dY_Over_dX : Data_Vector    := (others => Zero);
  
     I_Start  : Index renames X_Stuff.I_Start;
     I_Finish : Index renames X_Stuff.I_Finish;
     S : X_Structure  renames X_Stuff;
     F : Coefficients renames Spline.F;
     V, Z   : Tri.Column;
  begin
     if not S.Initialized then
        raise Must_Call_Prepare_X_Data_First;
     end if;
 
     -- The Spline variable must carry info about the start and
     -- and end of the calculated spline.  So this stuff was input
     -- in the call to Prepare_X_Data, placed into variable
     -- X_Stuff. Now it's placed in variable Spline, so that it
     -- it can be passed on to the interpolation procedures.  Its
     -- complicated, but less error prone this way.
 
     Spline.I_Start  := I_Start;
     Spline.I_Finish := I_Finish;
 
     for i in I_Start .. I_Finish-1 loop
        dY(i) := Y_Data(i+1) - Y_Data(i);
     end loop;
 
     for i in I_Start .. I_Finish-1 loop
        dY_Over_dX(i)  := dY(i) * S.dX_Inverse(i);
     end loop;
 
     -- Solve the matrix equation M*Z = V for vector Z.  First get V:
     --
     -- The following boundary conditions differ in the cases of
     -- Natural Splines, and Clamped Splines.
     --
     -- set Boundary_Value_First = Alpha*Y_dot_dot + Beta*Y_dot to get:
     --
     --  (2*h_0*Beta - 6*Alpha)*C_0 + Beta*h_0*C_1
     --         = -3*Boundary_Value_First + 3*Beta*(Y_1 - Y_0)/h_0.
     --
     -- At the other end set Boundary_Value_Last = Alpha2*Y_dot_dot + Beta2*Y_dot,
     -- to get:
     --
     --  (2*h_n-1*Beta2 + 6*Alpha2)*C_n + Beta2*h_n-1*C_n-1
     --         = 3*Boundary_Value_Last - 3*Beta2*(Y_n - Y_n-1)/h_n-1.
 
     V(I_Start)   := Three * (S.Bound_First.Beta * dY_Over_dX(I_Start)
                          - S.Bound_First.Boundary_Val);
 
     V(I_Finish) := Three * (-S.Bound_Last.Beta * dY_Over_dX(I_Finish-1)
                          + S.Bound_Last.Boundary_Val);
 
     -- Get rest of matrix M and vector V.
     -- Both clamped and Natural_BC Spline types have the following in common.
  
     for i in I_Start+1 .. I_Finish-1 loop
        V(i) := Three * (dY_Over_dX(i) - dY_Over_dX(i-1));
     end loop;
 
     -- Solve for Z in the equation M*Z = V:
     -- Matrix M has been overwritten with the LU decomposition of M.
 
     Tri.Solve (Z, S.M, V, I_Start, I_Finish);
 
     -- coefficients of X**2:
     for i in I_Start .. I_Finish loop
        F(2)(i) := Z(i);
     end loop;
 
     -- coefficients of X**3:
     F(3)(I_Finish) := Zero;
     for i in I_Start .. I_Finish-1 loop
        F(3)(i) := (F(2)(i+1) - F(2)(i)) * S.dX_Inverse(i) * One_Third;
     end loop;
 
     -- coefficients of X**1:
     for i in I_Start .. I_Finish-1 loop
        F(1)(i) := - S.dX(i) * (F(2)(i+1) + Two * F(2)(i)) * One_Third;
     end loop;
 
     for i in I_Start .. I_Finish-1 loop
        F(1)(i) := F(1)(i) + dY_Over_dX(i);
     end loop;
 
     -- Following is used only for testing:
     F(1)(I_Finish) := F(1)(I_Finish-1) +
                     S.dX(I_Finish-1)*(F(2)(I_Finish-1) + 
                     F(2)(I_Finish));
 
     -- F : Coefficients renames Spline.F;
     -- So Spline.F is now initialized
 
  end Get_Spline;
 
  --------------
  -- Value_At --
  --------------
 
  function Value_At 
    (X      : in Real;
     X_Data : in Data_Vector;
     Y_Data : in Data_Vector;
     Spline : in Spline_Coefficients) 
     return Real
  is
     S, Segment : Index;
     Y, dX1, dX2, dX3 : Real;
 
     F : Coefficients renames Spline.F;
     I_Start  : Index renames Spline.I_Start;
     I_Finish : Index renames Spline.I_Finish;
  begin
 
     if X > X_Data(I_Finish) or X < X_Data(I_Start) then  
        raise Constraint_Error with "Independent variable X out of range";
     end if;
 
     -- Find the spline segment associated with independent variable X.
     -- "Segment" is the Index of data array X_Data such that
     --    X_Data(Segment) <= X < X_Data(Segment+1).
 
     for Test_Segment in I_Start+1 .. I_Finish loop
        if X < X_Data(Test_Segment) then
           Segment := Test_Segment-1;
           exit;
        end if;
     end loop;
 
     if X = X_Data(I_Finish) then
        Segment := I_Finish-1;
     end if;
 
     -- sum the polynomial, return Y
 
     S   := Segment;
     dX1 := X - X_Data(Segment);
     dX2 := dX1 * dX1;
     dX3 := dX1 * dX2;
 
     Y := Y_Data(S) + F(1)(S)*dX1 + F(2)(S)*dX2 + F(3)(S)*dX3;
 
     return Y;
 
  end Value_At;
 
  -------------------------
  -- First_Derivative_At --
  -------------------------
 
  function First_Derivative_At 
    (X      : in Real;
     X_Data : in Data_Vector;
     Spline : in Spline_Coefficients) 
     return Real
  is
     S, Segment : Index;
     dY_over_dX, dX1, dX2 : Real;
 
     F : Coefficients renames Spline.F;
     I_Start  : Index renames Spline.I_Start;
     I_Finish : Index renames Spline.I_Finish;
  begin
 
     if X > X_Data(I_Finish) or X < X_Data(I_Start) then  
        raise Constraint_Error with "Independent variable X out of range";
     end if;
 
     -- Find the spline segment associated with independent variable X.
     -- "Segment" is the Index of data array X_Data such that
     --    X_Data(Segment) <= X < X_Data(Segment+1).
 
     for Test_Segment in I_Start+1 .. I_Finish loop
        if X < X_Data(Test_Segment) then
           Segment := Test_Segment-1;
           exit;
        end if;
     end loop;
 
     if X = X_Data(I_Finish) then
        Segment := I_Finish-1;
     end if;
 
     -- Integrate each polynomial segment from a1=X_Data(i) to b1=X_Data(I+1).
     -- y = y_1 + F_1 * (x - a1) + F_2 * (x - a1)**2 + F_3 * (x - a1)**3
 
     -- sum the polynomial, return dY / dX
 
     S   := Segment;
     dX1 := X - X_Data(Segment);
     dX2 := dX1 * dX1;
 
     dY_over_dX := F(1)(S) + Two * F(2)(S)*dX1 + Three * F(3)(S)*dX2;
 
     return dY_over_dX;
 
  end First_Derivative_At;
 
  --------------------------
  -- Second_Derivative_At --
  --------------------------
 
  function Second_Derivative_At 
    (X      : in Real;
     X_Data : in Data_Vector;
     Spline : in Spline_Coefficients) 
     return Real
  is
     S, Segment : Index;
     ddY, dX1 : Real;
 
     F : Coefficients renames Spline.F;
     I_Start  : Index renames Spline.I_Start;
     I_Finish : Index renames Spline.I_Finish;
  begin
     -- Independent_Variable_X_Out_Of_Range:
 
     if X > X_Data(I_Finish) or X < X_Data(I_Start) then  
        raise Constraint_Error;
     end if;
 
     -- Find the spline segment associated with independent variable X.
     -- "Segment" is the Index of data array X_Data such that
     --    X_Data(Segment) <= X < X_Data(Segment+1).
 
     for Test_Segment in I_Start+1 .. I_Finish loop
        if X < X_Data(Test_Segment) then
           Segment := Test_Segment-1;
           exit;
        end if;
     end loop;
 
     if X = X_Data(I_Finish) then
        Segment := I_Finish-1;
     end if;
 
     -- Integrate each polynomial segment from a1=X_Data(i) to b1=X_Data(I+1).
     -- y = y_1 + F_1 * (x - a1) + F_2 * (x - a1)**2 + F_3 * (x - a1)**3
 
     -- sum the polynomial, return dY / dX
 
     S   := Segment;
     dX1 := X - X_Data(Segment);
 
     ddY := Two * F(2)(S) + Six * F(3)(S)*dX1;
 
     return ddY;
 
  end Second_Derivative_At;
 
  --------------
  -- Integral --
  --------------
 
  function Integral 
    (X_Data : in Data_Vector;
     Y_Data : in Data_Vector;
     Spline : in Spline_Coefficients) 
     return Real
  is
     Area, Segment_Area : Real;
     a1, b1 : Real;
     d1, d2, d3, d4 : Real;
 
     F : Coefficients renames Spline.F;
     I_Start  : Index renames Spline.I_Start;
     I_Finish : Index renames Spline.I_Finish;
     Half    : constant Real := +0.5;
     Quarter : constant Real := +0.25;
  begin
     -- Integrate each polynomial segment from a1=X_Data(i) to b1=X_Data(I+1).
     -- y = y_1 + F_1 * (x - a1) + F_2 * (x - a1)**2 + F_3 * (x - a1)**3
 
     Area := Zero;
 
     for i in I_Start .. I_Finish-1 loop
        a1 := X_Data(i);
        b1 := X_Data(i+1);
        d1 := b1 - a1;
        d2 := d1 * d1;
        d3 := d1 * d2;
        d4 := d2 * d2;
 
        Segment_Area :=  
           Y_data(i) * d1
           + F(1)(i) * d2 * Half
           + F(2)(i) * d3 * One_Third
           + F(3)(i) * d4 * Quarter; 
 
        Area := Area + Segment_Area;
     end loop;
 
     return Area;
 
  end Integral;
 
end Spline;

