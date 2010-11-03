/**
   \page page6 Profiling MAOS 

   MAOS is a pthread enabled software, and in the current stage, the gprof based
   profiling does not generate useful information. Instead, we have to reply on
   system wide profilers, such as the oprofile.

   Prepare:

   <code>opcontrol --callgraph=16</code><p>
   <code>opcontrol --start </code><p>
   
   Then Run MAOS
   
   View profile:

   <code>opcontrol --stop</code><p>
   <code>opcontrol --dump</code><p>
   
   <code>opreport -cgf -l maos | gprof2dot.py -f oprofile | dot -Tpng -o output.png</code><p>

   or 

   <code>opgprof maos | gprof2dot.py | dot -Tpng -o output.png</code><p>

   gprof2dot.py can be found in the scripts folder.
*/
