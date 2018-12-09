
-- Spline_tst2 tests spline calculations on a subset of the Data_Vector.

with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;
with Spline;

procedure Spline_tst_2 is

   type Real is digits 15;

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use Math;
   package flio is new Float_IO (Real);
   use flio;

   type Index is range 1  ..  128;

   type Data_Vector is array(Index) of Real;

   package Sin_Spline is new Spline (Real, Index, Data_Vector);
   use Sin_Spline;

   X_Data, Y_Data : Data_Vector := (others => 0.0);
   S : Spline_Coefficients;

   -- tst2 tests spline calculations on a subset of the Data_Vector:

   I_start  : constant Index := Index'First + 2;
   I_finish : constant Index := Index'Last - 2;

   Pii    : constant Real := Ada.Numerics.Pi;
   DeltaX : Real := Pii / (Real (I_finish) - Real (I_start));
   X, Y   : Real;
   Y_true  : Real;
   dY, ddY : Real;
   DeltaX1, DeltaX2, DeltaX3 : Real;
   Integral_of_Sin, Err : Real;

   Bound_First, Bound_Last : Boundary;
   X_Stuff : X_Structure;

   -- Requires Text_IO.

   procedure Pause (s1,s2,s3,s4,s5,s6 : string := " ") is
      Continuation : Character := ' ';
   begin
      Text_IO.New_Line;
      if S1 /= " " then put_line (S1); end if;
      if S2 /= " " then put_line (S2); end if;
      if S3 /= " " then put_line (S3); end if;
      if S4 /= " " then put_line (S4); end if;
      if S5 /= " " then put_line (S5); end if;
      if S6 /= " " then put_line (S6); end if;
      Text_IO.New_Line;
      begin
         Text_IO.Put ("Type a character to continue: ");
         Text_IO.Get_Immediate (Continuation);
      exception
         when others => null;
      end;
      Text_IO.New_Line;
   end pause;

begin

   --  Start by making natural splines

   DeltaX := Pii / (Real (I_finish) - Real (I_start));

   for I in Index loop
      X         := DeltaX * (Real (I) - Real(Index'First));
      X_Data(I) := X;
      Y_Data(I) := Sin (X);
   end loop;

   --  Natural boundary conditions:

   Bound_First := (Alpha => 1.0, Beta => 0.0, Boundary_Val => 0.0);
   Bound_Last  := (Alpha => 1.0, Beta => 0.0, Boundary_Val => 0.0);

   Prepare_X_Data (X_Stuff, X_Data, I_start, I_finish,
                   Bound_First, Bound_Last);

   Get_Spline (S, X_Stuff, Y_Data);

   --  Check that Spline is continuous

   Pause
    ("Test 1: Check that the Spline is continuous at the knots.",
     "Recall that the spline equations were derived via this assumption",
     "(along with the assumption that the first 2 derivatives of the curve",
     "are also continuous at the knots.)");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      Y := Y_Data(I) + S.F(1)(I)*DeltaX1 + S.F(2)(I)*DeltaX2 + S.F(3)(I)*DeltaX3;
      -- This is the Y at I+1 due to spline starting at point I.
      -- Compare with Y at I+1 = Y_Data(I+1)

      Put (Y - Y_Data(I+1)); New_Line;

   end loop;

   --  Check that derivative is continuous

   Pause ("Verify that the derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      dY := S.F(1)(I) + 2.0*S.F(2)(I)*DeltaX1 + 3.0*S.F(3)(I)*DeltaX2;
      -- This is the dY at I+1 due to spline starting at point I.
      -- Compare with dY at I+1 = S.F(1)(I+1)

      Put (dY - S.F(1)(I+1)); New_Line;

   end loop;

   --  Check that 2nd derivative is continuous

   Pause ("Verify that the 2nd derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      ddY := 2.0*S.F(2)(I) + 6.0*S.F(3)(I)*DeltaX1;
      -- This is the ddY at I+1 due to spline starting at point I.
      -- Compare with ddY at I+1 = 2.0*S.F(2)(I+1)

      Put (ddY - 2.0*S.F(2)(I+1)); New_Line;

   end loop;

   Pause ("Subtract Spline prediction from true value.");

   for I in I_start .. I_finish-1 loop
      X     := X_Data(I) + DeltaX/2.0;
      Y     := Value_At (X, X_Data, Y_Data, S);
      Y_true := Sin (X);
      Put (Y - Y_true); New_Line;
   end loop;
   New_line;

   --  Check quadrature:

   Pause
    ("Test 1b: Check numerical integration of area under the curve.");

   Integral_of_Sin := -(Cos (X_Data(I_finish)) - Cos (X_Data(I_start)));
   Err := Abs (Integral_of_Sin - Integral (X_Data, Y_Data, S));

   New_Line;
   Put ("Error in numerical quadrature ="); Put (Err); 
 --Put (Integral_of_Sin); Put (Integral (X_Data, Y_Data, S));
   New_Line;

   --  make natural splines, test Cos.

   DeltaX := Pii / (Real (I_finish) - Real (I_start));

   for I in I_start .. I_finish loop
      X         := DeltaX * (Real (I) - Real(I_start));
      X_Data(I) := X;
      Y_Data(I) := Cos (X);
   end loop;

   --  Natural boundary conditions:
   Bound_First := (Alpha => 1.0, Beta => 0.0, Boundary_Val => -1.0);
   Bound_Last  := (Alpha => 1.0, Beta => 0.0, Boundary_Val => 1.0);

   Prepare_X_Data (X_Stuff, X_Data, I_start, I_finish,
                   Bound_First, Bound_Last);

   Get_Spline (S, X_Stuff, Y_Data);

   --Check that Spline is continuous
   Pause
    ("Test 1c: Check that the Spline is continuous at the knots.",
     "In this test the boundary conditions are given by values of",
     "the second derivative at the end points, but this second derivative",
     "is not 0.0 as in the case of true natural splines.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      Y := Y_Data(I) + S.F(1)(I)*DeltaX1 + S.F(2)(I)*DeltaX2 + S.F(3)(I)*DeltaX3;
      -- This is the Y at I+1 due to spline starting at point I.
      -- Compare with Y at I+1 = Y_Data(I+1)

      Put (Y - Y_Data(I+1)); New_Line;

   end loop;

   --  Check that derivative is continuous

   Pause ("Verify that the derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      dY := S.F(1)(I) + 2.0*S.F(2)(I)*DeltaX1 + 3.0*S.F(3)(I)*DeltaX2;
      -- This is the dY at I+1 due to spline starting at point I.
      -- Compare with dY at I+1 = S.F(1)(I+1)

      Put (dY - S.F(1)(I+1)); New_Line;

   end loop;

   --  Check that 2nd derivative is continuous

   Pause ("Verify that the 2nd derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      ddY := 2.0*S.F(2)(I) + 6.0*S.F(3)(I)*DeltaX1;
      -- This is the ddY at I+1 due to spline starting at point I.
      -- Compare with ddY at I+1 = 2.0*S.F(2)(I+1)

      Put (ddY - 2.0*S.F(2)(I+1)); New_Line;

   end loop;

   Pause ("Subtract Spline prediction from true value.");

   for I in I_start .. I_finish-1 loop
      X     := X_Data(I) + DeltaX/2.0;
      Y     := Value_At (X, X_Data, Y_Data, S);
      Y_true := Cos (X);
      Put (Y - Y_true); New_Line;
   end loop;
   New_line;

   -- Make mixed splines: natural at one end, clamped at the other.

   DeltaX := 1.5 * Pii / (Real (I_finish) - Real (I_start));

   for I in I_start .. I_finish loop
      X         := DeltaX * (Real (I) - Real(I_start));
      X_Data(I) := X;
      Y_Data(I) := Sin (X);
   end loop;

   -- Natural:
   Bound_First := (Alpha => 1.0, Beta => 0.0, Boundary_Val => 0.0);

   -- Clamped with derivative = 0.0:
   Bound_Last  := (Alpha => 0.0, Beta => 1.0, Boundary_Val => 0.0);

   Prepare_X_Data (X_Stuff, X_Data, I_start, I_finish,
                     Bound_First, Bound_Last);

   Get_Spline (S, X_Stuff, Y_Data);

   -- Check that Spline is continuous:

   Pause
    ("Test 2: Check that the mixed-boundary Spline is continuous at the knots.",
     "Recall that the spline equations were derived via this assumption",
     "(along with the assumption that the first 2 derivatives of the curve",
     "are also continuous at the knots.)");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      Y := Y_Data(I) + S.F(1)(I)*DeltaX1 + S.F(2)(I)*DeltaX2 + S.F(3)(I)*DeltaX3;
      -- This is the Y at I+1 due to spline starting at point I.
      -- Compare with Y at I+1 = Y_Data(I+1)

      Put (Y - Y_Data(I+1)); New_Line;

   end loop;

   --  Check that derivative is continuous
   Pause ("Make sure that the derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      dY := S.F(1)(I) + 2.0*S.F(2)(I)*DeltaX1 + 3.0*S.F(3)(I)*DeltaX2;
      -- This is the dY at I+1 due to spline starting at point I.
      -- Compare with dY at I+1 = S.F(1)(I+1)

      Put (dY - S.F(1)(I+1)); New_Line;

   end loop;

   --  Check that 2nd derivative is continuous
   Pause ("Make sure that the 2nd derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      ddY := 2.0*S.F(2)(I) + 6.0*S.F(3)(I)*DeltaX1;
      -- This is the ddY at I+1 due to spline starting at point I.
      -- Compare with ddY at I+1 = 2.0*S.F(2)(I+1)

      Put (ddY - 2.0*S.F(2)(I+1)); New_Line;

   end loop;


   Pause ("Subtract Spline prediction from true value.");

   for I in I_start .. I_finish-1 loop
      X     := X_Data(I) + DeltaX/2.0;
      Y     := Value_At (X, X_Data, Y_Data, S);
      Y_true := Sin (X);
      Put (Y - Y_true); New_Line;
   end loop;
   New_line;


   -- Make mixed splines: natural at one end, clamped at the other.
   -- This time we get them wrong.


   DeltaX := 1.5 * Pii / (Real (I_finish) - Real (I_start));

   for I in I_start .. I_finish loop
      X        := DeltaX * (Real (I) - Real(I_start));
      X_Data(I) := X;
      Y_Data(I) := Sin (X);
   end loop;

   -- Now we have the boundary conditions wrong for the Sin curve:
   -- Natural:
   Bound_First := (Alpha => 0.0, Beta => 1.0, Boundary_Val => 0.0);

   -- Clamped with derivative = 1.0:
   Bound_Last  := (Alpha => 1.0, Beta => 0.0, Boundary_Val => 0.0);

   Prepare_X_Data (X_Stuff, X_Data, I_start, I_finish,
                   Bound_First, Bound_Last);

   Get_Spline (S, X_Stuff, Y_Data);

  --Check that Spline is continuous
  Pause
   ("Test 3: Make sure the mixed-boundary Spline is continuous at the knots.",
    "Recall that the spline equations were derived via this assumption",
    "(along with the assumption that the first 2 derivatives of the curve",
    "are also continuous at the knots.)");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      Y := Y_Data(I) + S.F(1)(I)*DeltaX1 + S.F(2)(I)*DeltaX2 + S.F(3)(I)*DeltaX3;
      -- This is the Y at I+1 due to spline starting at point I.
      -- Compare with Y at I+1 = Y_Data(I+1)

      Put (Y - Y_Data(I+1)); New_Line;

   end loop;

   --  Check that derivative is continuous
   Pause ("Make sure that the derivative of the spline is continuous.");

      for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      dY := S.F(1)(I) + 2.0*S.F(2)(I)*DeltaX1 + 3.0*S.F(3)(I)*DeltaX2;
      -- This is the dY at I+1 due to spline starting at point I.
      -- Compare with dY at I+1 = S.F(1)(I+1)

      Put (dY - S.F(1)(I+1)); New_Line;

   end loop;

   --  Check that 2nd derivative is continuous
   Pause ("Make sure that the 2nd derivative of the spline is continuous.");

   for I in I_start .. I_finish-1 loop

      DeltaX1 := X_Data(I+1) - X_Data(I);
      DeltaX2 := DeltaX1*DeltaX1;
      DeltaX3 := DeltaX1*DeltaX2;

      ddY := 2.0*S.F(2)(I) + 6.0*S.F(3)(I)*DeltaX1;
      -- This is the ddY at I+1 due to spline starting at point I.
      -- Compare with ddY at I+1 = 2.0*S.F(2)(I+1)

      Put (ddY - 2.0*S.F(2)(I+1)); New_Line;

   end loop;

   Pause
     ("Subtract Spline prediction from true value.",
      "In this test we deliberately made the Boundary conditions",
      "wrong, so the answers should be bad at the end points.");

   for I in I_start .. I_finish-1 loop
      X      := X_Data(I) + DeltaX/2.0;
      Y      := Value_At (X, X_Data, Y_Data, S);
      Y_true := Sin (X);
      Put (Y - Y_true); New_Line;
   end loop;

   New_line;

end;
