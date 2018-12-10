
-- package Spline
--
-- Package of cubic splines, for interpolation, integration,
-- and differentiation of tabulated data. 
--
-- Splines of this sort are especially useful for smooth data 
-- sets (Y_1, Y_2, ... Y_n) defined at unequally spaced 
-- knots: X_1, X_2, X_3, ... X_n. They are not least squares
-- fits, so they are not much use for fitting noisy data.
-- 
-- The splines are, optionally, Natural or Clamped, at either end 
-- point. The user inputs a value for either the derivative of the 
-- curve of the second derivative of the curve at the end points. If
-- the user sets the second derivative of the curve to 0.0 at the
-- end point, then the spline (at that end point) is called Natural.
-- (But the value of the 2nd derivative there can be anything.) If
-- instead he chooses to input a value for the 1st derivative at the
-- end point then the spline (at that end point) is called clamped.
--
-- These spline routines are for interpolation, not curve
-- fitting, so they are not useful for data with noise. (The
-- interpolated curve produced passes through all of the data
-- points input.) The user inputs points (X1,Y1) (X2,Y2)..(Xn,Yn).
-- The X points are input as X_Data(i), the Y points as Y_Data(i).
-- The spline interpolation algorithm creates a function F(X) such
-- that (X, F(X)) passes through each the points input.
-- The function F(X) can then be used predict the value of Y
-- associated with a value of X that lies between Xm and Xm+1.
--
-- The interpolated curve is a different 3rd order polynomial for
-- each Xm.  It gives the value of Y at point X near Xm according to:
--
--   Y = F(X) = Ym + F(1,m)*(X-Xm) + F(2,m)*(X-Xm)^2 + F(3,m)*(X-Xm)^3,
--
-- assuming that X satisfies Xm <= X < Xm+1.  The Spline coefficients
-- F(O,m) are calculated by solving equations
-- derived from the requirement that the spline segments (given above)
-- are continuous at each Xm, and have the same first and second
-- derivatives at each Xm, (with the exception of m = n).  However,
-- at the first X and the last X, the spline segments may have
-- zero valued second derivatives (Natural Splines).  The user also
-- has the option of inputting the values for the first derivatives
-- at the end points (Clamped Splines).
--
-- Notes on algorithm
--
-- A cubic spline is a collection of third order polynomials S(X)
-- connecting tabulated data points (X_i, Y_i):
--
--   S_i(X) = Y_i + B_i * (X-X_i) + C_i * (X-X_i)**2 + D_i * (X-X_i)**3
--
-- Here the knots X_i are indexed on the range 0, 1, ..., N.  So the
-- spline segments are valid on the range i = 0,..N-1.
-- In what follows, X_i+1 - X_i will be called h_i, the standard variable
-- in numerical analysis texts for DeltaX.  (In the code below these
-- deltas will be called dX, and (B, C, D) will be (F(1), F(2), F(3)).)
-- From the requirement that the polynomials S_i are continuous at the
-- knots X = X_i, and have continuous 1st and 2nd derivatives there, you
-- get the simultaneous equations,
--
--    Y_i+1 = Y_i + B_i * h_i + C_i * h_i**2 + D_i * h_i**3   (1)
--
--    B_i+1 = B_i + 2 * C_i * h_i + 3 * D_i * h_i**2     (2)
--
--    C_i+1 = C_i + 3 * D_i * h_i         (3)
--
-- These are valid on the range i = 0,..,N-1.  At the end points,
-- X_0 and X_N, there are no equations; we must impose some value on
-- the 1st or the 2nd derivative of Y there.
-- So these are the equations to solve for B, C, and D, given Y and
-- h_i.  To solve these we eliminate B and D from the equations, and
-- solve for C.  (We also need to imput boundary conditions into the
-- the equations, by fixing either the 1st of 2nd derivative of Y at
-- the end points.)  First eliminate D from equations (1) and (2) by
-- solving for D in (3), and substituting into (1) and (2):
--
--    Y_i+1 = Y_i + B_i * h_i + (2*C_i +C_i+1) * h_i**2 / 3    (4)
--
--    B_i+1 = B_i + (C_i + C_i+1) * h_i        (5)
--
-- Again, i = 0,..,N-1.
-- Finally, solve for B_i in (4), reduce the Indices by one, and plug
-- into (5) to get equations solely in terms of C_i:
--
-- h_j-1*C_j-1 + 2(h(j-1 + h_j)*C_j + h_j*C_j+1 = 3*dY_i/h_i - 3*dY_i-1/h_i-1
--
-- where dY_i = Y_i+1 - Y_i and where i is in 1,...,N-1.  So we now
-- have N-1 equations for N+1 variables C_i, where i is in 0,...,N.
-- The last two equations will come from boundary conditions at
-- X_0 and X_N.  Remember that C is twice the second derivative of Y
-- so that if we impose a value on the second derivative of Y at the
-- end points, (call it Y_dot_dot), then we have the final two equations:
-- C_0 = Y_dot_dot_First / 2, and C_N = Y_dot_dot_Last / 2.  Suppose
-- instead we wish to impose a value on the 1st derivative of Y at
-- one or both end points.  (Call it Y_dot.)  This is trickier. The
-- equations are B_0 = Y_dot_First, and B_N = Y_dot_Last.  The full
-- equations are in terms of the variable C, so we must eliminate B from
-- the above two with equ. (4) for the i=0 end point and equ. (5)
-- substituted into (4) for the i=n end point.  At i=0 we get,
--
--  Y_dot_First = B_0 = (Y_i+1 - Y_i)/h_i - (2*C_i + C_i+1)*h_i/3  (i=0)
--
--  Y_dot_Last  = B_n = B_n-1 + (C_n-1 + C_n)*h_n-1  (i=n)
--
-- The first equation above is
--
--   2*h_0*C_0 + h_0*C_1 = -3*Y_dot_First + 3*(Y_1 - Y_0)/h_0
--
-- The second equation above is plugged into (4) at i=n-1 to give
--
-- h_n-1*C_n-1 + 2*h_n-1*C_n = 3*Y_dot_Last - 3*(Y_n - Y_n-1)/h_n-1.
--
-- So at either end point you can specify Y_dot or Y_dot_dot to get
-- unique solutions of the equations that establish continuity of the
-- cubic polynomial spline segments and their first two derivatives.
--
-- Sometimes mixed boundary conditions are required.  Instead of specifying
-- Y_dot or Y_dot_dot, you impose a value on Alpha*Y_dot_dot + Beta*Y_dot.
-- (A common reason for doing this is that the solutions of some PDE
-- equations have boundary conditions of this sort.)  It should
-- be clear how to do this now.  The two alternatives at the first end
-- point are:
--
--   2*h_0*C_0 + h_0*C_1 = -3*Y_dot_First + 3*(Y_1 - Y_0)/h_0,
--
--          C_0 = Y_dot_dot_first / 2.
--
-- Multiply the second equation by -6*Alpha, multiply the first
-- equation by Beta, add them together to get:
--
--  (2*h_0*Beta - 6*Alpha)*C_0 + Beta*h_0*C_1
--         = -3*(Alpha*Y_dot_dot_First + Beta*Y_dot_First) +
--           + 3*Beta*(Y_1 - Y_0)/h_0.
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
--
-- (Might not always find solutions for arbitrary values of Alpha, Beta,
-- and Boundary_Value.)
--

with Tridiagonal_LU;

generic

  type Real is digits <>;
  type Index is range <>;
  type Data_Vector is array(Index) of Real;

package Spline is

   type Coefficients is array(1..3) of Data_Vector;

   type Spline_Coefficients is record
      F        : Coefficients;
      I_Start  : Index;
      I_Finish : Index;
   end record;

   -- For each boundary point, the following boundary condition is defined:
   --
   --     Alpha * Y_dot_dot + Beta * Y_dot = Boundary_Val.
   --
   -- If Alpha = 1 and Beta = 0 then the boundary condition is that the
   -- the 2nd derivative of the curve Y_dot_dot is set to value
   -- Boundary_Val.  (if then Boundary_Val is 0.0, then this is called a
   -- Natural spline.)  If Alpha = 0 and Beta = 1 then the boundary condition
   -- is that the 1st derivative of the curve Y_dot is set to a value
   -- Boundary_Val. This is called a clamped spline.  In the two special
   -- cases given above, a unique spline exists and can be calculated by
   -- the routines below.  In some cases mixed boundary conditions are
   -- required, Alpha and Beta both non-zero.  In this case we can't
   -- guarantee that a unique spline can be found satisfying this boundary
   -- condition and satisfying the equations of continuity.  Use with
   -- care under these circumstances.  This option is provided because
   -- some differential equations satisfy mixed boundary conditions.

   type Boundary is record
      Alpha : Real        := 1.0;
      Beta  : Real        := 0.0;
      Boundary_Val : Real := 0.0;
   end record;

   Natural_BC : constant Boundary :=
     (Alpha        => 1.0, 
      Beta         => 0.0,
      Boundary_Val => 0.0);

   -- Procedure Prepare_X_Data 
   -- prepares the data Arrays for use by Get_Spline. 
   -- Its a separate procedure so that it can be removed
   -- from inner loops .. might call this once, and Get_Spline many 
   -- times. The procedure only cares about the X positions of the
   -- knots (X(i)). In many problems the knots remain constant, but
   -- the Y values (Y(i)) change in an inner loop.
   -- The procedure also performs the LU decomposition in preparation
   -- for solution of the coupled equations that determine the
   -- spline coefficients.

   type X_Structure is limited private;

   procedure Prepare_X_Data 
     (X_Stuff     :    out X_Structure;
      X_Data      : in     Data_Vector;
      I_Start     : in     Index    := Index'first;
      I_Finish    : in     Index    := Index'Last;
      Bound_First : in     Boundary := Natural_BC;
      Bound_Last  : in     Boundary := Natural_BC);

   -- Procedure Get_Spline 
   -- calculates the coefficients of powers of X in the cubic
   -- spline polynomial: F(1), F(2), and F(3). (These F's are
   -- in the record Spline.)  If Ym and and Xm are elements of 
   -- the set of data points being fit with a spline (the knots), 
   -- then a value of the function Y where    Ym <= Y < Ym+1 
   -- is given by the interpolation formula:
   --
   --  Y = Ym  +  F(1,m)*(X-Xm)  +  F(2,m)*(X-Xm)^2  +  F(3,m)*(X-Xm)^3,
   --
   -- Procedure Get_Spline only calculates Spline.F.  To get Y at points
   -- not equal to Ym, Ym=1 etc, use function Value_At.
   --
   -- Notice that the 1st derivatives of the spline at the data points
   -- (Xm, Ym) are given by F(1,m), the second derivatives are given
   -- by 2*F(2,m) and the 3rd derivative by 6*F(3,m).  To get derivatives
   -- of the spline away from the data points (knots), one must use
   -- the interpolatory formulas derived from the equation for Y above.

   procedure Get_Spline 
     (Spline  :    out Spline_Coefficients;
      X_Stuff : in     X_Structure;
      Y_Data  : in     Data_Vector);

   function Value_At 
     (X      : in Real;
      X_Data : in Data_Vector;
      Y_Data : in Data_Vector;
      Spline : in Spline_Coefficients) 
      return Real;

   function First_Derivative_At 
     (X      : in Real;
      X_Data : in Data_Vector;
      Spline : in Spline_Coefficients) 
      return Real;

   function Second_Derivative_At 
     (X      : in Real;
      X_Data : in Data_Vector;
      Spline : in Spline_Coefficients) 
      return Real;
   -- Second_Derivative_At is highly inaccurate.

   function Integral 
     (X_Data : in Data_Vector;
      Y_Data : in Data_Vector;
      Spline : in Spline_Coefficients) 
      return Real;

   Must_Call_Prepare_X_Data_First : exception;

private

   package Tri is new Tridiagonal_LU (Real, Index);

   type X_Structure is record
      Initialized : Boolean     := False;
      dX          : Data_Vector := (others => 0.0);
      dX_Inverse  : Data_Vector := (others => 0.0);
      M           : Tri.Matrix  := (others => (others => 0.0));
      I_Start     : Index       := Index'First;
      I_Finish    : Index       := Index'Last;
      Bound_First : Boundary    := Natural_BC;
      Bound_Last  : Boundary    := Natural_BC;
   end record;

end Spline;
