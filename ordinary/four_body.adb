
with Ada.Numerics.Generic_Elementary_Functions;

package body Four_Body is

   One : constant Real := +1.0;
   Two : constant Real := +2.0;
   Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin + 32);

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use mth;

   procedure Update_State
     (Y          : in out Dynamical_Variable;
      Body_id    : in     Bodies;
      X, Z, U, W : in     Real)
   is
      State_id : State_Index;
   begin
      State_id := (Body_id - Bodies'First) * No_Of_Bodies;
      Y(State_id+0) := X;
      Y(State_id+1) := Z;
      Y(State_id+2) := U;
      Y(State_id+3) := W;
   end Update_State;

   function State_Val
     (Y       : Dynamical_Variable;
      Body_id : Bodies;
      XYUV_id : XYUV_Index)
      return Real
   is
      State_id : constant State_Index 
              := (Body_id - Bodies'First) * No_Of_Bodies + XYUV_id;
   begin
      return Y(State_id);
   end State_Val;

   --  Define a metric for the vector Y

   function Norm (Y : Dynamical_Variable)
      return Real
   is
      Sum : Real := +0.0;
      X, Z : Real;
   begin
      for Body_id in Bodies loop
         X   := State_Val (Y, Body_id, 0);
         Z   := State_Val (Y, Body_id, 1);
         Sum := Sum + Abs (X) + Abs (Z);
      end loop;
      return Sum;
   end Norm;

   function "-"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable
   is
      Result : Dynamical_Variable;
   begin
      for I in State_Index loop
         Result(I) := Left(I) - Right(I);
      end loop;
      return Result;
   end "-";

   function "*"
     (Left  : Real;
      Right : Dynamical_Variable)
      return Dynamical_Variable
   is
      Result : Dynamical_Variable;
   begin
     for X in State_Index loop
        Result(X) := Left * Right(X);
     end loop;
     return Result;
   end "*";

   function "+"
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable)
      return Dynamical_Variable
   is
      Result : Dynamical_Variable;
   begin
      for X in State_Index loop
         Result(X) := Left(X) + Right(X);
      end loop;
     return Result;
   end "+";

   --  Defines the differential equation to integrate:
  --
   --   dY/dt = F (t, Y)
   --
   --  The N-Body equation of motion is:
   --
   --  Acceleration(Body_j)  = (2 * Pi)**2 *
   --           SumOver(Body_k) { Mass(Body_k) * DeltaR / NORM(DeltaR)**3 },
   --
   --  where DeltaR = (X(Body_k) - X(Body_j)),
   --  and mass, time, and distance are in the natural units given above.
   --  Actually, the (2 * Pi)**2 is absorbed into the Mass to avoid
   --  unecessary multiplications.
   --
   --  An optimization: each mass in a pair A and B experiences a force that
   --  is equal in magnitude to the force the other mass experiences in the
   --  pair, so there's no need to calculate the forces twice.

   function F
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable
   is
      Deriv : Dynamical_Variable;
      Delta_X, Delta_Z  : Real;
      Rinv, R2, R3, MR3 : Real;
      X2, X0, Z2, Z0    : Real;
      State_id, Other_State_id : State_Index;
   begin
      for The_Body in Bodies loop
         State_id := (The_Body - Bodies'First) * No_Of_Bodies;
         Deriv(State_id + 0) := Y(State_id + 2);
         Deriv(State_id + 1) := Y(State_id + 3);
         Deriv(State_id + 2) := 0.0;
         Deriv(State_id + 3) := 0.0;
      end loop;

      for The_Body in Bodies loop

         State_id := (The_Body - Bodies'First) * No_Of_Bodies;

         for Other_Body in Bodies loop
         if Other_Body > The_Body then

            --  The Other_Bodies with indices < The_Body
            --  already had their contributions
            --  to the force on The_Body added to the sum.

            Other_State_id := (Other_Body - Bodies'First) * No_Of_Bodies;

            X2     := State_Val (Y, Other_Body, 0);
            X0     := State_Val (Y, The_Body, 0);
            Z2     := State_Val (Y, Other_Body, 1);
            Z0     := State_Val (Y, The_Body, 1);
            Delta_X  := X2 - X0;
            Delta_Z  := Z2 - Z0;
         -- Delta_X  :=  Y(Other_Body)(0) - Y(The_Body)(0);
         -- Delta_Z  :=  Y(Other_Body)(1) - Y(The_Body)(1);
            R2       :=  Delta_X*Delta_X + Delta_Z*Delta_Z + Min_Allowed_Real;
            Rinv     :=  Sqrt (One / R2);
            R3       :=  Rinv * Rinv * Rinv;

            --  force on The_Body due to the other bodies with ID > The_Body
            MR3 :=  Mass (Other_Body) * R3;
            Deriv(State_id + 2) :=  Deriv(State_id + 2) + Delta_X*MR3;
            Deriv(State_id + 3) :=  Deriv(State_id + 3) + Delta_Z*MR3;

            --  force on Other_Body due to the The_Body. (Accel. = force/Mass)
            MR3 :=  Mass (The_Body) * R3;
            Deriv(Other_State_id + 2) :=  Deriv(Other_State_id + 2) - Delta_X*MR3;
            Deriv(Other_State_id + 3) :=  Deriv(Other_State_id + 3) - Delta_Z*MR3;

         end if;
         end loop;
      end loop;

      return Deriv;

   end F;

end Four_Body;
