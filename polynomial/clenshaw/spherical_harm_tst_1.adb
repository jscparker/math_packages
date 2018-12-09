
with Spherical_Harmonics;
with gauss_quadrature_61;
with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;

procedure Spherical_Harm_Tst_1 is

   type Real is digits 15;
   
   package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real); 
   use Maths;

   type Base_Poly_Index is new Integer;
   
   subtype Poly_Index is Base_Poly_Index range 0..1000;
   
   package Sph is new Spherical_Harmonics 
     (Real, Sqrt, Exp, Log, Base_Poly_Index, Poly_Index'Last);
   use Sph;
   
   package rio is new Float_IO(Real);
   use rio;

   package Quad is new Gauss_Quadrature_61 (Real, Sqrt);
   use Quad;
   subtype Gauss_Index is Quad.Gauss_Index_61;


   Num_Points : constant := 32 * 16; 
   -- need 32 * 256 to do l=3000  (m=small)
   -- need 32 * 8   to do l=700   (m=small)

   type X_axis is range 1 .. Num_Points;

   DeltaX : constant Real := (X_Upper_Bound - X_Lower_Bound) / Real (Num_Points);

   Area, Max, Diff, Y, W, Y1, Y2 : Real;
   
   X_Start, X_Final : Real;
   Error, d_Area    : Real;
   F_val : Function_Values;
   X_g   : Gauss_Values;
   
   l_max, m : Real := 1.0;
   L, L1, L2, Int_L_Max : Poly_Index;
   Int_m : Poly_Index;
   k_max : Poly_Index;

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8 : string := "") is
      Continue : Character := ' ';
   begin
      new_line;
      if S1 /= "" then put_line (S1); end if;
      if S2 /= "" then put_line (S2); end if;
      if S3 /= "" then put_line (S3); end if;
      if S4 /= "" then put_line (S4); end if;
      if S5 /= "" then put_line (S5); end if;
      if S6 /= "" then put_line (S6); end if;
      if S7 /= "" then put_line (S7); end if;
      if S8 /= "" then put_line (S8); end if;

      dialog: loop
         begin
            New_Line;
            Put ("Type a character to continue: ");
            Get_Immediate (Continue);
            exit dialog;
         exception
            when others => null;
         end;
      end loop dialog;
      new_line;
   end pause;

begin
   
   --
   -- Test 1. Compare the two calculations of Polys
   --
   New_Line;
   Pause
    ("Test 1: Compare two independent calculations of Spherical Harmonics.",
     "The difference between Spherical_Harm_2 and Spherical_Harm is printed.",
     "Spherical_Harm_2 is used only for testing, and can fail if m>144, so",
     "test values of m should be limited to less than 144.");
         
   new_line;
   put_line ("Enter the value of the magnetic quantum number m (0, 1, ..), e.g. 40: ");
   get (m);
   new_line;
   put ("Enter the max value of the azimuthal quantum number l (m, m+1, ..), e.g. 60: ");
   get (l_max);

   Int_l_max := Poly_Index (l_max);
   Int_m     := Poly_Index (m);
   k_max     := Int_l_max - Int_m;

   for k in Poly_Index range 0 .. k_max loop

      L   := k + Int_m;
      Max := 0.0;

      for I in X_Axis loop
         X_start := X_Lower_Bound + DeltaX * (Real (I) - Real (X_Axis'First));
         X_Final := X_Start + DeltaX;
         --  Interval of Gaussian quadrature.

         Find_Gauss_Nodes (X_Start, X_Final, X_g);
         for N in Gauss_Index loop
            Y1   := Spherical_Harm_2 (L, m, X_g(N));   -- use for testing, m>144 fails.
            Y2   := Spherical_Harm (L, m, X_g(N));     -- works ok for high L, m
            Diff := Abs (Y1 - Y2) / (Abs(Y1) + Abs (Y2) + 1.0e-16);
            if (Diff > Max) then
               Max := Diff;
            end if;
         end loop;
      end loop;

      new_line; 
      put ("L = "); put (Poly_Index'Image(L));
      put(" Max difference between the two = "); put (Max);

   end loop;
   
   --
   -- Test 2. Create polynomials, check Norm
   -- k goes from 0...k_max, and l automatically goes from m, m+1...m+k...
   --
   New_Line;
   Pause 
    ("Test 2: test normalization of the polynomials.",
     "The polynomials will be squared and integrated on (-1,1),",
     "and the result, will be printed in column 2 below.");

   new_line;
   put_line ("Enter the value of the magnetic quantum number m (0, 1, ..): ");
   get (m);
   new_line;
   put ("Enter the max value of the azimuthal quantum number l (m, m+1, ..): ");
   get (l_max);

   Int_l_max := Poly_Index (l_max);
   Int_m     := Poly_Index (m);
   k_max     := Int_l_max - Int_m;

   for k in Poly_Index range 0 .. k_max loop

      L    := k + Int_m;
      Area := 0.0;

      for I in X_Axis loop
         X_start := X_Lower_Bound + DeltaX * (Real (I) - Real (X_Axis'First));
         X_Final := X_Start + DeltaX;
         --  Interval of Gaussian quadrature.

         Find_Gauss_Nodes (X_Start, X_Final, X_g);
         for N in Gauss_Index loop
            W         := Poly_Weight (X_g(N));
            Y         := Spherical_Harm (L, m, X_g(N));
            F_val (N) := W*Y*Y;
         end loop;
         Get_Integral (F_val, X_Start, X_Final, d_Area, Error);
         Area := Area + d_Area;
      end loop;

      new_line; put ("L = "); put (Poly_Index'Image(L));
      put(" Norm = "); put (Area);

   end loop;

   New_Line;
   Pause 
    ("Test 3: test orthogonality of the Spherical Harmonics.",
     "Two polynomials will be multiplied and integrated on (-1,1)",
     "and the result, will be printed in column 2.",
     "One of the Polynomials will be fixed at l = m.",
     "The other polynomials will be in the range l = m..l_max.");
         
   new_line;
   put_line ("Enter the value of the magnetic quantum number m (0, 1, ..): ");
   get (m);
   new_line;
   put ("Enter l_max (m, m+1, ..): ");
   get (l_max);
   new_line;

   Int_l_max := Poly_Index (l_max);
   Int_m     := Poly_Index (m);
   k_max     := Int_l_max - Int_m;

   for k in Poly_Index range 0 .. k_max loop

      L1   := Int_m;
      L2   := k + Int_m;
      Area := 0.0;

      for I in X_Axis loop
         X_start := X_Lower_Bound + DeltaX * (Real (I) - Real (X_Axis'First));
         X_Final := X_Start + DeltaX;
         --  Interval of Gaussian quadrature.

         Find_Gauss_Nodes (X_Start, X_Final, X_g);
         for N in Gauss_Index loop
            W         := Poly_Weight (X_g(N)); -- Check this
            Y1        := Spherical_Harm (L1, m, X_g(N));
            Y2        := Spherical_Harm (L2, m, X_g(N));
            F_val (N) := W*Y1*Y2;
         end loop;
         Get_Integral (F_val, X_Start, X_Final, d_Area, Error);
         Area := Area + d_Area;
      end loop;

      new_line;
      put ("L2 = "); put (Poly_Index'Image(L2));
      put("  Inner product of the 2 polys = "); put (Area);

   end loop;

end Spherical_Harm_Tst_1;
