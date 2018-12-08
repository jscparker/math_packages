
--  Demonstrates use of package:   Disorderly.Random

with Disorderly.Random; use Disorderly.Random; 
with Disorderly.Random.Clock_Entropy;
with Text_io; use text_io;

procedure random_demo_1 is

   X : Random_Int;

   Stream_1 : State;  
   --  Must declare one of these for each desired independent stream of rands.
   --  To get successive random nums X from stream_k, 
   --  you then call Get_Random(X, stream_k);

begin

   new_line;
   put_line ("The generator needs an initial state before it can create streams");
   put_line ("of random numbers.  We usually call the state Stream_1.  You can");
   put_line ("create lots of independent Streams by calling Reset with different");
   put_line ("values of the seeds, which are called Initiator1, Initiator2 ... Below");
   put_line ("is a series of 12 states created by procedure Reset. The 12 states are");
   put_line ("created by setting Initiator1 to 1, 2, 3 ... 12, while keeping ");
   put_line ("the other 3 Initiators constant. In fact we only needed to change");
   put_line ("Initiator1 by just 1 bit to get a complete change in all 4 of the");
   put_line ("64-bit Integers comprising the state. Each line below shows the four 64 bit");
   put_line ("integers of a state. Twelve independent states are printed on 12 lines.");
   declare Continue : Character;
   begin
      new_line; put  ("Enter a character to continue: ");
      get_immediate (Continue);
   exception
      when others => null;
   end;

   for k in Seed_Random_Int range 1..12 loop
     Reset (Stream_1, k, 333, 4444, 55555);
     new_line; put (Image (Stream_1));
   end loop;

   new_line(2);
   put_line ("Next we test functions Image and Value. Function Image translates");
   put_line ("the State (an array of 4 Integers) into a String. Function Value");
   put_line ("translates the string back to array of 4 Integers.");
   put_line ("We do this back and forth, and print the results below. Each string");
   put_line ("of 4 numbers (representing a State) should appear twice.");
   declare Continue : Character;
   begin
      new_line; put  ("Enter a character to continue: ");
      get_immediate (Continue);
   exception
      when others => null;
   end;

   for k in Seed_Random_Int range 1..4 loop
     Reset (Stream_1, k, 1, 1);
     new_line; put (Image (Stream_1));
     new_line; put (Image (Value (Image (Stream_1))));
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
   declare Continue : Character;
   begin
      new_line; put  ("Enter a character to continue: ");
      get_immediate (Continue);
   exception
      when others => null;
   end;

   -- up top we wrote: use Disorderly.Random;
   -- so we can just call Clock_Entropy.Reset instead of
   -- Disorderly.Random.Clock_Entropy.Reset:

   new_line;
   for k in Seed_Random_Int range 1..12 loop
   --Disorderly.Random.Clock_Entropy.Reset (Stream_1);
     Clock_Entropy.Reset (Stream_1);
     new_line; put (Image (Stream_1));
   end loop;


   new_line(2);
   put_line ("Finally, print 8 random nums from Stream_1, the state just initialized.");
   declare Continue : Character;
   begin
      new_line; put  ("Enter a character to continue: ");
      get_immediate (Continue);
   exception
      when others => null;
   end;

   new_line;
   for k in 1..8 loop
     Get_Random (X, Stream_1);
     new_line; put (Random_Int'Image (X));
   end loop;

end;

