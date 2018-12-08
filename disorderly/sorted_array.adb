
--with Text_io; use Text_io;
package body Sorted_Array is

   Empty_Slot : constant Item := Item'Last; 
   --  make type Item bigger if possible than data *set*.
   --  eg Item = Parent_Random_Int; or Random_Int'Base. 
   --  would prefer Empty_Slot is out of range of actual data.
   --  not necessary tho.

   type Table is array (Table_Index) of Item;

   T : Table := (others => Empty_Slot);

   Index_of_Next_Item_to_Read : Table_Index := Table_Index'First;

   No_More_Items_Remaining : Boolean := False;

   --  Routines for reading the array sequentially (from beginning to end):

   --------------------------------------
   -- Start_Reading_Array_at_Beginning --
   --------------------------------------

   procedure Start_Reading_Array_at_Beginning is
   begin
     Index_of_Next_Item_to_Read := Table_Index'First;
   end Start_Reading_Array_at_Beginning;

   -------------------------------
   -- Set_Index_to_Present_Item --
   -------------------------------

   --  procedure Set_Index_to_Present_Item leaves index alone if it points
   --  a readable item, or else goes to next item if it exists, or to end of array
   --  if there are no more items.  So if after calling Set_Index_to_Next_Item we
   --  find that T(Index_of_Next_Item_to_Read) = Empty_Slot then that means we have
   --  read every item in array already: set No_More_Items_Remaining := True.

   procedure Set_Index_to_Present_Item is
   begin
      if T(Index_of_Next_Item_to_Read) /= Empty_Slot then
         return;
         -- Index is set to proper next Item. Leave it here.
      end if;
      for i in table_Index loop
         if Index_of_Next_Item_to_Read < Table_Index'Last then
            Index_of_Next_Item_to_Read := Index_of_Next_Item_to_Read + 1;
	 else
           No_More_Items_Remaining := True;
	   return; -- Index_of_Next_Item_to_Read = Table_Index'Last, and no Item found
	 end if;
	 if T(Index_of_Next_Item_to_Read) /= Empty_Slot then
           --  Index is set to proper next Item. Leave it here.
	   return;
	 end if;
      end loop;
   end Set_Index_to_Present_Item;

   procedure Get_Next_Item
     (X : out Item;
      Item_is_Invalid : out Boolean)
   is
   begin
      -- init out param: 

      Item_is_Invalid := False;

      Set_Index_to_Present_Item; 
      -- stay at present index unless empty slot; skip all possible empty slots

      if No_More_Items_Remaining then -- go home early.
         Item_is_Invalid :=  True;
         X := Empty_Slot;
	 return;
      end if;

      X := T(Index_of_Next_Item_to_Read);

      --  if we got to Table_Index'Last then Item at Table_Index'Last was valid, 
      --  but now there may be no Items left:

      if Index_of_Next_Item_to_Read = Table_Index'Last then
         No_More_Items_Remaining := True; -- the reason we need No_More_Items_Re..
      else
         --  Set_Index so Next Item is found next time around:
         Index_of_Next_Item_to_Read := Index_of_Next_Item_to_Read + 1;
      end if;

   end Get_Next_Item;




   -- Overflow Stacks:

   subtype Stack_Index is Table_Index range 0 .. Size_of_Overflow_Stacks - 1;
   type Item_Array is array (Stack_Index) of Item;

   type Stack is
   record
     Storage : Item_Array := (others => Empty_Slot);
     Next_Free_Slot_id : Stack_Index := Stack_Index'First;
   end record;

   Collision_Stack : Stack;


   type Sorted_Stack is
   record
     Storage : Item_Array := (others => Empty_Slot);
     Next_Free_Slot_id : Stack_Index := Stack_Index'First;
   end record;

   Low_End_Stack, High_End_Stack : Sorted_Stack;



   function No_of_Items_in_Low_Stack return Table_Index is
   begin
     return Table_Index (Low_End_Stack.Next_Free_Slot_id);
   end No_of_Items_in_Low_Stack;

   -----------------------------------
   -- Initialize_Table_for_Restart --
   -----------------------------------

   procedure Initialize_Table_for_Restart is
   begin

      T := (others => Empty_Slot);
      Collision_Stack.Storage  := (others => Empty_Slot);
      Low_End_Stack.Storage    := (others => Empty_Slot);
      High_End_Stack.Storage   := (others => Empty_Slot);
      Collision_Stack.Next_Free_Slot_id := Stack_Index'First;
      Low_End_Stack.Next_Free_Slot_id   := Stack_Index'First;
      High_End_Stack.Next_Free_Slot_id  := Stack_Index'First;
      No_More_Items_Remaining           := False;
      Index_of_Next_Item_to_Read        := Table_Index'First;

   end Initialize_Table_for_Restart;

   ----------
   -- Push --
   ----------

   procedure Push
     (X        : in Item;
      Stack_id : in out Stack)
   is 
      pragma Assert (Stack_Index'First = 0);
      --  So No of items in stack = Next_Free_Slot_id
      S : Item_Array renames Stack_id.Storage;
      Next_Free_Slot_id : Stack_Index renames Stack_id.Next_Free_Slot_id;
   begin

      S (Next_Free_Slot_id) := X;
 
      if Next_Free_Slot_id = Stack_Index'Last then raise Constraint_Error; end if;

      Next_Free_Slot_id := Next_Free_Slot_id + 1;

   end Push;
   

   -------------------
   -- Push_and_Sort --
   -------------------

   -- only way to update a Sorted_Stack, or sort fails.

   procedure Push_and_Sort
     (X        : in Item;
      Stack_id : in out Sorted_Stack)
   is
      pragma Assert (Stack_Index'First = 0);
      X_Item_Size : constant Item := X; 
      tmp : Item; 
      S : Item_Array renames Stack_id.Storage;
      Next_Free_Slot_id : Stack_Index renames Stack_id.Next_Free_Slot_id;
   begin

      S (Next_Free_Slot_id) := X;
 
      bubble_sort:
      for i in reverse Stack_Index'First+1 .. Next_Free_Slot_id loop

        if X_Item_Size < S (i-1) then  -- X is now in S(Next_Free_Slot_id)
	   tmp     := S (i-1);
	   S (i-1) := S (i);
	   S (i)   := tmp;
        elsif X_Item_Size = S (i-1) then
	   Push (X, Collision_Stack);
	   exit bubble_sort;
	else -- X_Item_Size > S (i-1)
	   exit bubble_sort;
	end if;

      end loop bubble_sort;

      if Next_Free_Slot_id = Stack_Index'Last then raise Constraint_Error; end if;

      Next_Free_Slot_id := Next_Free_Slot_id + 1;

   end Push_and_Sort;

   ---------------------
   -- Insert_and_Sort --
   ---------------------

   -- X is both the item to be inserted, and the Item_Size

   procedure Insert_and_Sort  
     (X : in Item)
   is
      X_Item_Size : constant Item := X; 
      X_index_0   : constant Table_Index  
                     := Table_Index (X_Item_Size / Items_per_Table_Entry);
      Empty_Slot_id : Table_Index;
      Empty_Slot_Detected : Boolean := False;

   begin

      if X = Empty_Slot then                  
         --  value Empty_Slot signifies empty space in array, so is reserved.
         Push_and_Sort (X, High_End_Stack);
         return;
      end if;

      if T (X_index_0) = Empty_Slot then
         T (X_index_0) := X;
         return;
      end if;

      if X_Item_Size = T(X_index_0) then

         Push (X, Collision_Stack);
         return;

      elsif X_Item_Size < T (X_index_0) then

        if X_index_0 = Table_Index'First then
           Push_and_Sort (X, Low_End_Stack);
           return;
	end if;

        Empty_Slot_Detected := False;
	Find_Empty_Slot_Lower_Down:
        for i in reverse Table_Index range Table_Index'First .. X_index_0-1 loop
	   if T(i) = Empty_Slot then
              Empty_Slot_id := i;
              Empty_Slot_Detected := True;
              exit Find_Empty_Slot_Lower_Down;
	   end if;
	end loop Find_Empty_Slot_Lower_Down;

        if Empty_Slot_Detected = False then
	   Push_and_Sort (X, Low_End_Stack);
	   return;
	end if;

        --  shift the empty slot back to the rt place, and fill it with X:

        -- common short cut:
	if Empty_Slot_id = X_index_0-1 then -- put X into the table
	   T(Empty_Slot_id) := X;
	   return;
	end if;

        Shift_Slot_Up:
        for i in Empty_Slot_id+1 .. X_index_0 loop -- i-1 is the empty slot.
           if X_Item_Size > T(i) then 
             T(i-1) := T(i);      -- shift T(i) into the Empty Slot (at i-1)
	   --T(i)   := Empty_Slot;-- will fill this below with X
	   elsif X_Item_Size = T(i) then
	     Push (X, Collision_Stack);
	     T(i-1) := X;
	     return;
	   else
	     T(i-1) := X;
	     return;
           end if;
	end loop Shift_Slot_Up;

      elsif X_Item_Size > T (X_index_0) then

        if X_index_0 = Table_Index'Last then
	   Push_and_Sort (X, High_End_Stack);
	   return;
	end if;

        Empty_Slot_Detected := False;
	Find_Empty_Slot_Higher_Up:
        for i in Table_Index range X_index_0+1 .. Table_Index'Last loop
	   if T(i) = Empty_Slot then
              Empty_Slot_id := i;
              Empty_Slot_Detected := True;
              exit Find_Empty_Slot_Higher_Up;
	   end if;
	end loop Find_Empty_Slot_Higher_Up;

        if Empty_Slot_Detected = False then
           Push_and_Sort (X, High_End_Stack);
           return;
	end if;


        --  shift the empty slot back to the rt place, and fill it with X:

        -- common short cut:
	if Empty_Slot_id = X_index_0+1 then -- put X into the table
	   T(Empty_Slot_id) := X;
	   return;
	end if;

        Shift_Slot_Down:
        for i in reverse X_index_0 .. Empty_Slot_id-1 loop -- i+1 is the empty slot.
           if X_Item_Size < T(i) then 
             T(i+1) := T(i);      -- put T(i) into the Empty Slot (at i+1)
	   --T(i)   := Empty_Slot;-- will fill this below with X 
	   elsif X_Item_Size = T(i) then
	     Push (X, Collision_Stack);
	     T(i+1) := X;
	     return;
	   else
	     T(i+1) := X;
	     return;
           end if;
	end loop Shift_Slot_Down;

      end if;

   end Insert_and_Sort;

   ---------------------------
   -- Array_Sort_Successful --
   ---------------------------

   function Array_Sort_Successful return Boolean is
      Hi, Lo : Item;
      Non_Empty_Slot_id : Table_Index;
   begin

      Find_Non_Empty_Slot_Higher_Up:
        for i in Table_Index loop
           if T(i) /= Empty_Slot then
              Non_Empty_Slot_id := i;
              Lo := T(i);
              exit Find_Non_Empty_Slot_Higher_Up;
	   end if;
	end loop Find_Non_Empty_Slot_Higher_Up;

        for i in Non_Empty_Slot_id+1 .. Table_Index'Last loop 
           if  T(i) = Empty_Slot then 
	     null;  
	   else
	     Hi := T(i);
	     if Hi < Lo then
	       return False;
	     end if;
	     Hi := Lo;
           end if;
	end loop;

        return True;

   end Array_Sort_Successful;

   -------------------------------
   -- No_of_Collisions_Detected --
   -------------------------------

   function No_of_Collisions_Detected return Table_Index is
   begin
     return Collision_Stack.Next_Free_Slot_id - Stack_Index'First; 
     --  subtype of Table_Index
   end No_of_Collisions_Detected;

end Sorted_Array;
