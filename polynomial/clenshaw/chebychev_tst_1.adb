
with Clenshaw;
with Chebychev; 
with gauss_quadrature_61;
with chebychev_quadrature;
with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;

-- Use data structure in Chebychev to instantiate Clenshaw, and make 
-- Chebychev polynomials of the 2nd kind.
-- Chebychev polynomials should be normalized respect a weight W(x).
-- In other words, Integral (Poly(x)*Poly(x)*W(x)) = 1.0
-- Test should show that Chebychev Quadrature works flawlessly on
-- Chebychev polynomials; Gaussian Quadrature less so; Summing
-- rectangles is worst of all.

procedure Chebychev_Tst_1 is

   type Real is digits 15;

   package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real); 
   use Maths;

   package Chebychev_Polys is new Chebychev (Real, Exp, Log);
   use Chebychev_Polys;

   Poly_Limit : constant Poly_ID_Integer := 200;

   package Cheby is new 
      Clenshaw (Real, Poly_ID_Integer, Poly_Limit, Alpha, Beta, Q_0);
   use Cheby;

   package realio is new Float_IO(Real);
   use Realio;

   package Gaussian_Quad is new gauss_quadrature_61 (Real);
   use Gaussian_Quad;

   package Cheby_Quad is new Chebychev_Quadrature (Real, 2**6, Cos, Sin);

   -- Divide [-1,1] on X axis into 2**n parts:

   No_of_Intervals : constant := 2**0; -- power of 2 best
   type X_axis is range 1 .. No_of_Intervals;

   DeltaX : constant Real := (X_Upper_Bound - X_Lower_Bound) / Real (No_of_Intervals);

   W, X, Y : Real;
   Q  : Poly_Values;

   d_Area, Area, Error : Real;
   X_Start, X_Final : Real;
   X_gauss : Gauss_Values;
   F_val   : Function_Values;
   -- gaussian quadrature over segments of (-1, 1).
   
   X_gauss_Cheby : Cheby_Quad.Gauss_Values;
   F_val_Cheby   : Cheby_Quad.Function_Values;

begin
   --  Sum rectangles for area:

   new_line(2); put ("Check normalization. Sum rectangles for area:");
   new_line(2); 
  
   Sum_Rectangles:
   declare
      No_of_Subdivisions : constant Positive := 256;
      dX : constant Real := DeltaX / Real (No_of_Subdivisions);
   begin 
      for k in Poly_Id_Type range Poly_Id_Type'First .. Poly_Id_Type'First+20 loop
   
         Area := 0.0;
         X    := X_Lower_Bound + 0.5 * dX;
   
         for m in 1 .. No_of_Subdivisions loop
         for i in X_Axis loop
            Evaluate_Qs (X, Q, k);
            Y     := Q(k);
            W     := Poly_Weight (X);
            Area  := Area + W*Y*Y;
            X     := X + dX;
         end loop;
         end loop;
   
         new_line; put ("k = "); put (Poly_Id_Type'Image(k));
         put("    Norm = "); put (Area * dX / Norm (k));
   
      end loop;
   end Sum_Rectangles;
   
   new_line(2); 
   put ("Check normalization. Use gaussian quadrature for area: "); 
   new_line(1); 
   
   for k in Poly_Id_Type range Poly_Id_Type'First .. Poly_Id_Type'First+30 loop

      Area := 0.0;

      for i in X_Axis loop
         X_start := X_Lower_Bound + DeltaX * (Real (i) - Real (X_Axis'First));
         X_Final := X_Start + DeltaX;
         --  Interval of Gaussian quadrature.

         Find_Gauss_Nodes (X_Start, X_Final, X_gauss);
         for N in Gauss_Index_61 loop
            Evaluate_Qs (X_gauss(N), Q, k);
            Y         := Q(k);
            W         := Poly_Weight (X_gauss(N));
            F_val (N) := W*Y*Y;
         end loop;
         Get_Integral (F_val, X_Start, X_Final, d_Area, Error);
         Area := Area + d_Area;
      end loop;

      new_line; put ("k = "); put (Poly_Id_Type'Image(k));
      put("    Norm = "); put (Area / Norm (k));

   end loop;
   
   new_line(2); 
   put ("Check normalization. Use chebychev quadrature for area: "); 
   new_line(1); 
   
   for k in Poly_Id_Type range Poly_Id_Type'First .. Poly_Id_Type'First+32 loop

      Area := 0.0;

      for i in X_Axis loop
         X_start := X_Lower_Bound + DeltaX * (Real (i) - Real (X_Axis'First));
         X_Final := X_Start + DeltaX;
         --  Interval of Gaussian quadrature.

         Cheby_Quad.Find_Gauss_Nodes (X_Start, X_Final, X_gauss_Cheby);
         for N in Cheby_Quad.Gauss_Index loop
            Evaluate_Qs (X_gauss_Cheby(N), Q, k);
            Y         := Q(k);
            W         := Poly_Weight (X_gauss_Cheby(N));
            F_val_Cheby (N) := W*Y*Y;
         end loop;
         Cheby_Quad.Get_Integral (F_val_Cheby, X_Start, X_Final, d_Area);
         Area := Area + d_Area;
      end loop;

      new_line; put ("k = "); put (Poly_Id_Type'Image(k));
      put("    Norm = "); put (Area / Norm (k));

   end loop;
   
end Chebychev_Tst_1;

