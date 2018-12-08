
--  At the moment the overflow stacks raise the constraint_error 
--  if the array is full..
--
--  Quickly written mess, but all that matters is that the output is
--  correctly sorted, so we just test this each time the routine is run.
--  
--  For a given X <= Max_Size_of_Item, we get its putative 
--  table index from X_Index := Size_of_X / S where S is some constant.
--  require      Max_Size_of_Item / S <=     Table_Index'Last.
--  or      1 +  Max_Size_of_Item / S <= 1 + Table_Index'Last.
--  using  Max_Size_of_Item <= S * (1 + Max_Size_of_Item / S):
--         Max_Size_of_Item <= S * (1 + Table_Index'Last).
--  Finally, S >= 1 + Max_Size_of_Item / (1 + Table_Index'Last).
--
--
-- put 10 items (the numbers 0..9) into 4 bins (with bin_id 0..3):
--
-- Item is an x below: 0..Max_Size_of_Item-1 = 0..9,    9=Max_Size_of_Item
--
--   |xxx|xxx|xxx|x00
--
-- Number of x's  =  10 = 1 + Max_Size_of_Item = Max_No_of_Items
-- Number of bins =  4  = 1 + Table_Index'Last = Size_of_Table
-- Items_per_bin  =  3  = 1 + (Max_No_of_Items-1) / No_of_Bins    = S
-- Items_per_bin  =  3  = 1 + (Max_No_of_Items-1) / Size_of_Table = S
--
-- bin_id = X / S = X / Items_per_bin


--  Table_Index'Last is supposed to be >> Max_Allowed_No_of_Items
--  Amazingly, equality of the 2 seems to work efficiently.  I don't
--  know why its efficient so I don't recommend it.
--
--  Typical case:
--
--   type Item is mod 2**63;
--
--   Max_Size_of_Item        := 2**48-1;
--   Max_Allowed_No_of_Items := 2**26;
--
--   type Table_Index is mod (2**26 + 2**24);

generic

   type Item is mod <>;
   --  can be 64 bit int, but max value that the array can store
   --  is Max_Size_of_Item.  (Here its 2**48-1).
   --  The Items are the Random ints produced by the generator.

   Max_Size_of_Item        : Item; -- Here Max_Size_of_Item < 2**48
   Max_Allowed_No_of_Items : Item;

   type Table_Index is mod <>; 
   --  make big enough to hold:  0..Max_Allowed_No_of_Items-1,
   --  but make it much bigger than that.

package Sorted_Array is

   pragma Assert (Table_Index'Last >= 
       Table_Index (Max_Allowed_No_of_Items) + Table_Index (Max_Allowed_No_of_Items/8));
   --  Table_Index'Last should be >> Max_Allowed_No_of_Items for best efficiency.
   --  Max_Allowed_No_of_Items/4  instead of /8 up there is better.

   --  2 routines for updating the table:

   procedure Insert_and_Sort
     (X : in Item);
 
   procedure Initialize_Table_for_Restart;


   --  2 routines for reading sequentially from the start:

   procedure Start_Reading_Array_at_Beginning;

   procedure Get_Next_Item
     (X : out Item;
      Item_is_Invalid : out Boolean);


   --  for reading the result:

   function No_of_Collisions_Detected return Table_Index;


   function Array_Sort_Successful return Boolean; -- for testing

   function No_of_Items_in_Low_Stack return Table_Index;

   Size_of_Overflow_Stacks : constant  := 2**12;

   Items_per_Table_Entry : constant Item  
                  := 1 + Max_Size_of_Item / (Item (Table_Index'Last) + 1);

end Sorted_Array;
