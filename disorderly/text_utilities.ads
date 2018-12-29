
-- Test handling for Random and Basic_Rand. Rewritten to pass gnatprove.

generic

   No_Of_Seeds : Positive;
   type Parent_Random_Int is mod <>;

package Text_Utilities 
   with Spark_Mode => On
is
   pragma Pure;

   --No_Of_Seeds : constant Positive := 3;
   --type Parent_Random_Int is mod 2**64;

   pragma Assert (Parent_Random_Int'Last = 2**64-1);

   Rand_Image_Width : constant := 20; -- 2^64-1 = 18_446_744_073_709_551_615
   Max_Image_Width  : constant Positive := Rand_Image_Width * No_Of_Seeds;

   subtype Random_Int_String_Index is Positive range 1 .. Rand_Image_Width;
   subtype Random_Int_String is String (Random_Int_String_Index);
 
   subtype Character_Digit is Character range '0' .. '9';

  -- All state values State.X(i) are less than 2^64.
  --
  -- All state values State.X(i) are in a range that is correctly evaluated by
  -- function 'Value'.
  --
  -- The calculation uses the modular arithmetic on modular type Parent_Random_Int. 
  -- If the image string encodes a number out-of-range of Parent_Random_Int,
  -- (i.e. > 2^64-1), then of course the function can't return that value.

   function Value (Random_Int_Image : in Random_Int_String) return Parent_Random_Int
   with Pre => 
     (for all i in Random_Int_String'Range => Random_Int_Image(i) in Character_Digit);

   function Image (X : in Parent_Random_Int) return Random_Int_String
   with Pre =>
      Parent_Random_Int'Last = 2**63 + (2**63-1) and
      Random_Int_String'Length = 20;
 
end Text_Utilities;

