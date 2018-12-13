
--  Demonstrates use of package:   Disorderly.Random

with Disorderly.Random; use Disorderly.Random; 
with Disorderly.Random.Clock_Entropy;
with Text_io; use text_io;

procedure Random_demo_1 is

   X : Random_Int;

   Stream_1 : State;  
   --  Must declare one of these for each desired independent stream of rands.
   --  To get successive random nums X from stream_k, 
   --  you then call Get_Random(X, stream_k);

   procedure Pause is
      Continue : Character;
   begin
      new_line; put  ("Enter a character to continue: ");
      get_immediate (Continue);
      new_line;
   exception
      when others => null;
   end Pause;

begin

   new_line;
   put_line ("The generator needs an initial state before it can create streams");
   put_line ("of random numbers.  We usually call the state Stream_1.  You can");
   put_line ("create lots of independent Streams by calling Reset with different");
   put_line ("values of the seeds, which are called Initiator1, Initiator2 ... Below");
   put_line ("is a series of 12 states created by procedure Reset. The 12 states are");
   put_line ("created by setting Initiator1 to 1, 2, 3 ... 12, while keeping ");
   put_line ("the other 2 Initiators constant. In fact we only needed to change");
   put_line ("Initiator1 by just 1 bit to get a complete change in all four of the");
   put_line ("64-bit Integers comprising the state. Each line below shows the four 64 bit");
   put_line ("integers of a state. Twelve independent states are printed on 12 lines.");
   Pause;

   for k in Seed_Random_Int range 1..12 loop
     Reset (Stream_1, k, 4444, 55555, 666666);
     new_line; put (Formatted_Image (Stream_1));
   end loop;

   new_line(2);
   put_line ("Test functions Image and Value.");
   put_line ("Function Image translates the State (an array of 4 Integers) into a String.");
   put_line ("Function Value translates the string back to array of 4 Integers.");
   put_line ("Do this back and forth, and print the results below. Each string");
   put_line ("of 4 numbers (representing a State) should appear twice.");
   Pause;

   for k in Seed_Random_Int range 1..4 loop
     Reset (Stream_1, k, 1, 1, 1);
     new_line; put (Formatted_Image (Stream_1));
     new_line; put (Formatted_Image (Value (Image (Stream_1))));
     new_line;
   end loop;

   new_line(2);
   put_line ("Test of procedure:    Clock_Entropy.Reset (Stream_1)");
   put_line ("Procedure Clock_Entropy.Reset calls Calendar in an attempt to create");
   put_line ("a new and unique initial state each time Clock_Entropy.Reset is called.");
   put_line ("Below we call Clock_Entropy.Reset 12 times, in order to get 12 unique");
   put_line ("initial states. The 12 initial states (strings of 4 numbers) are printed");
   put_line ("on the 12 lines below. If you get the same State twice in a row then the");
   put_line ("procedure failed to find a new and unique initial state.");
   Pause;

   -- up top we wrote: use Disorderly.Random;
   -- so we can just call Clock_Entropy.Reset instead of
   -- Disorderly.Random.Clock_Entropy.Reset:

   new_line;
   for k in Seed_Random_Int range 1..12 loop
   --Disorderly.Random.Clock_Entropy.Reset (Stream_1);
     Clock_Entropy.Reset (Stream_1);
     new_line; put (Formatted_Image (Stream_1));
   end loop;


   new_line(2);
   put_line ("Print 8 random nums from Stream_1, the state just initialized.");
   Pause;

   new_line;
   for k in 1..8 loop
     Get_Random (X, Stream_1);
     new_line; put (Random_Int'Image (X));
   end loop;

   new_line(2);
   put_line ("Final Test");
   put_line ("Translate state Integers to a string, and then back to Integer. Do this");
   put_line ("back and forth for a long time. If error detected, then report failure.");
   Pause;

   Clock_Entropy.Reset (Stream_1);

   for k in Seed_Random_Int range 1 .. 2**26 loop
     Get_Random (X, Stream_1);
     if not Are_Equal (Stream_1, Value (Image (Stream_1))) then
        raise Program_Error with "FAILURE in Image/Value routines";
     end if;
   end loop;

   new_line; put  ("Finished part 1 of the final test.");

   for k in Seed_Random_Int range 1 .. 2**22 loop
     Reset (Stream_1, k, 1, 1, 1);
     if not Are_Equal (Stream_1, Value (Image (Stream_1))) then
        raise Program_Error with "FAILURE in Image/Value routines";
     end if;
   end loop;

   new_line; put  ("Finished part 2 of the final test.");

   for k in Seed_Random_Int range 1 .. 2**20 loop
     Clock_Entropy.Reset (Stream_1);
     if not Are_Equal (Stream_1, Value (Image (Stream_1))) then
        raise Program_Error with "FAILURE in Image/Value routines";
     end if;
   end loop;

   new_line; put  ("Finished part 3 of the final test.");

end Random_Demo_1;

