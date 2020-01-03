/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \page page60 Profiling MAOS 

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
