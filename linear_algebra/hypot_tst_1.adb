
with Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;
with Hypot;
with text_io; use text_io;

-- Test estimates err in calculation of Sqrt(a^2 + b^2), a/Sqrt(a^2 + b^2), etc.
--
-- Best accuracy uses  -mfpmath=387 ! Runs about the same speed at -mfpmath=387 and
-- -mfpmath=sse.

procedure hypot_tst_1 is

  type Real is digits 15;

  package Hypotenuse is new Hypot (Real); --use Hypotenuse;

  type e is digits 18;
  package Hypotenuse_e is new Hypot (e);

  package mathr is new Ada.Numerics.Generic_Elementary_Functions (Real);  use mathr;


  procedure Get_Hypotenuse_e
    (a, b           : in e;
     Hypot          : out e; 
     sn_e, cn_m_1_e : out e) renames Hypotenuse_e.Get_Hypotenuse;

  procedure Get_Hypotenuse_r
    (a, b : in Real;
     Hypot_er  : out e; 
     sn_er     : out e;
     cn_m_1_er : out e)
  is
     Hypot, cs_m_1, sn : Real;
  begin
     Hypotenuse.Get_Hypotenuse
       (a, b, Hypot, sn, cs_m_1);
     Hypot_er  := e (Hypot);
     sn_er     := e (sn);
     cn_m_1_er := e (cs_m_1);
  end Get_Hypotenuse_r;

  
  Zero : constant Real := +0.0;
  Two  : constant Real := +2.0;

  Min_Real : constant e := e (Two ** (Real'Machine_Emin + 2));


  a, b, hypot_0, sn_0, cn_m_1 : Real := Zero;

  Max_Err_0, Err_0, Max_Err_1, Err_1 : Real := Zero;
  Max_Err_2, Err_2, Max_Err_3, Err_3 : Real := Zero;
  Max_Err_4, Err_4, Max_Err_5, Err_5 : Real := Zero;
  Ave_Err_0, Ave_Err_1, Ave_Err_2, Ave_Err_3, Ave_Err_4, Ave_Err_5 : Real := Zero;
  
  sn_e, Hypot_e, cn_m_1_e : e;
  sn_er, Hypot_er, cn_m_1_er : e;

  Epsilon : constant Real := 1.00012345e-9 * 2.0;

  type int_64 is range 0 .. 2**63-1;
  No_Of_Bins         : constant int_64 := 1024;
  No_Of_Test_Vectors : constant int_64 := 2 * 10**8 / 4;

begin

   declare
      hypot_1, x_1, y_1 : Real := Zero;
   begin
      Hypotenuse.Get_Hypotenuse (0.0, 0.0, Hypot_1, x_1, y_1);
      put (Real'Image (Hypot_1));
      Hypotenuse.Get_Hypotenuse (1.0E-307, 1.0E-307, Hypot_1, x_1, y_1);
      put (Real'Image (Hypot_1));
      Hypotenuse.Get_Hypotenuse (1.0E-307, 1.0E+307, Hypot_1, x_1, y_1);
      put (Real'Image (Hypot_1));
      Hypotenuse.Get_Hypotenuse (1.0E-306, 1.0E-300, Hypot_1, x_1, y_1);
      put (Real'Image (Hypot_1));
   end;
   
   b := 0.000000000001 + 0.0;
   a := 0.881111111111;
   
   for j in 1 .. No_Of_Bins loop
   
      Max_Err_0 := Zero;
      Max_Err_1 := Zero;
      Max_Err_2 := Zero;
      Max_Err_3 := Zero;
      Max_Err_4 := Zero;
      Max_Err_5 := Zero;
   
      Ave_Err_0 := Zero;
      Ave_Err_1 := Zero;
      Ave_Err_2 := Zero;
      Ave_Err_3 := Zero;
      Ave_Err_4 := Zero;
      Ave_Err_5 := Zero;
   
      for i in int_64 range 1 .. No_Of_Test_Vectors loop
   
         Get_Hypotenuse_e
           (e(a), e(b), Hypot_e, sn_e, cn_m_1_e);
   
         -- uncorrected hypot:
   
         hypot_0 := Sqrt (a**2 + b**2);
       --hypot_e := Sqrt (e(a)**2 + e(b)**2);
   
         Err_0 := Real (Abs(e(Hypot_0) - Hypot_e) + Min_Real) 
              / (Abs(Hypot_0) + 1.0e22 * Real(Min_Real));
   
         if Err_0 > Max_Err_0 then Max_Err_0 := Err_0; end if;
         Ave_Err_0 := Ave_Err_0 + Err_0;
   
         -- Hypot_er is slightly corrected by Get_Hypotenuse_r:
   
         Get_Hypotenuse_r 
           (a, b, Hypot_er, sn_er, cn_m_1_er);

         -- test func Hypotenuse(a,b)
--         declare
--            hypot_1, tst_1, x_1, y_1 : Real := Zero;
--         begin
--            Hypotenuse.Get_Hypotenuse
--              (a, b, Hypot_1, x_1, y_1);
--            tst_1 := Hypotenuse.Hypotenuse (a, b);
--          --if Abs (tst_1 - Hypot_1) > 6.0E-17 then
--            if Abs (tst_1 - Hypot_1) > 0.0 then
--               put(real'image(tst_1 - Hypot_1));
--            end if;
--         end;
            
         Err_1 := Real (Abs(Hypot_er - Hypot_e) + Min_Real) 
              / (Abs(Hypot_0) + 1.0e22 * Real(Min_Real));
   
         if Err_1 > Max_Err_1 then Max_Err_1 := Err_1; end if;
         Ave_Err_1 := Ave_Err_1 + Err_1;
   
         -- uncorrected sn:
   
         sn_0 := Real'Min(a, b) / (hypot_0 + 2.0 * Real (Min_Real));
   
         Err_2 := Real (Abs(e(sn_0) - sn_e) + Min_Real) 
              / Real (Abs(sn_e) + 1.0e22 * Min_Real);
   
         if Err_2 > Max_Err_2 then Max_Err_2 := Err_2; end if;
         Ave_Err_2 := Ave_Err_2 + Err_2;
   
         -- sn_er is slightly corrected by Get_Hypotenuse_r:
   
         Err_3 := Real (Abs(sn_er - sn_e) + Min_Real) 
              / Real (Abs(sn_e) + 1.0e22 * Min_Real);
   
         if Err_3 > Max_Err_3 then Max_Err_3 := Err_3; end if;
         Ave_Err_3 := Ave_Err_3 + Err_3;
   
         -- uncorrected cn_m_1:
   
         cn_m_1 := - Real'Min(a, b)**2 / (hypot_0*(Real'Max(a, b) + hypot_0));
   
         Err_4 := Real (Abs(e(cn_m_1) - cn_m_1_e) + Min_Real) 
              / Real (Abs(cn_m_1_e) + 1.0e22 * Min_Real);
   
         if Err_4 > Max_Err_4 then Max_Err_4 := Err_4; end if;
         Ave_Err_4 := Ave_Err_4 + Err_4;
   
         -- cn_m_1_er is slightly corrected by Get_Hypotenuse_r:
   
         Err_5 := Real (Abs(cn_m_1_er - cn_m_1_e) + Min_Real) 
              / Real (Abs(cn_m_1_e) + 1.0e22 * Min_Real);
   
         if Err_5 > Max_Err_5 then Max_Err_5 := Err_5; end if;
         Ave_Err_5 := Ave_Err_5 + Err_5;
   
         b := b + Epsilon;
   
      end loop;
   
      --Epsilon := Epsilon + 2.7777e-10;
   
      new_line;
      put ("b =  "); put (Real'Image (b));
   
      new_line(2);
      put("Max Err in uncorrected hypot: "); put (Real'Image (Max_Err_0));
      new_line;
      put("Max Err in corrected hypot:   "); put (Real'Image (Max_Err_1));
      new_line(2);
      put("Max Err in uncorrected sin:   "); put (Real'Image (Max_Err_2));
      new_line;
      put("Max Err in sin:               "); put (Real'Image (Max_Err_3));
      new_line(2);
      put("Max Err in uncorrected cos-1: "); put (Real'Image (Max_Err_4));
      new_line;
      put("Max Err in cos-1:             "); put (Real'Image (Max_Err_5));
      new_line(2);
   
      new_line;
      put("Average Err in uncorrected hypot: "); 
      put (Real'Image (Ave_Err_0 / Real (No_Of_Test_Vectors)));
      new_line;
      put("Average Err in corrected hypot:   ");
      put (Real'Image (Ave_Err_1 / Real (No_Of_Test_Vectors)));
      new_line(2);
      put("Average Err in uncorrected sin:   ");
      put (Real'Image (Ave_Err_2 / Real (No_Of_Test_Vectors)));
      new_line;
      put("Average Err in sin:               ");
      put (Real'Image (Ave_Err_3 / Real (No_Of_Test_Vectors)));
      new_line(2);
      put("Average Err in uncorrected cos-1: ");
      put (Real'Image (Ave_Err_4 / Real (No_Of_Test_Vectors)));
      new_line;
      put("Average Err in cos-1:             ");
      put (Real'Image (Ave_Err_5 / Real (No_Of_Test_Vectors)));
      new_line(2);
   
   end loop;
   
end hypot_tst_1;
