
generic

   type Parent_Random_Int is mod <>;

package Xor_Shift
   with Spark_Mode => On
is 

   pragma Pure (Xor_Shift);
 
   pragma Assert (Parent_Random_Int'Modulus >= 2**61);
 
   subtype Valid_XOR_Shift_Range is Parent_Random_Int range 1 .. 2**61-1;
   -- Seeds need to be in this range. Input of 0 gives period of 1.
 
   procedure Get_Random_XOR_Shift_61 (Random_x : in out Parent_Random_Int);
   pragma Inline (Get_Random_XOR_Shift_61);
   -- Rearrangement generator. Takes numbers in the range 1 .. 2**61-1,
   -- and returns them in a different order. Has a period of 2**61-1, so
   -- each number in the range 1 .. 2**61-1 is returned exactly once during
   -- the 2**61-1 period. 
   --
   -- Random_x = 0 gives period of 1; needs to be rejected as a seed.
 
   procedure Get_Random_XOR_Shift_61_b (Random_x : in out Parent_Random_Int);
 
end Xor_Shift;
