/*
@ @licstart  The following is the entire license notice for the
JavaScript code in this file.

Copyright (C) 1997-2017 by Dimitri van Heesch

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

@licend  The above is the entire license notice
for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "MAOS", "index.html", [
    [ "Overview", "index.html", null ],
    [ "Introduction", "page10_intro.html", null ],
    [ "Compile the Code", "page20_compile.html", [
      [ "Compiling", "page20_compile.html#autotoc_md0", [
        [ "Prerequisites", "page20_compile.html#autotoc_md1", null ],
        [ "Download the code", "page20_compile.html#autotoc_md2", [
          [ "Option 1 (preferred)", "page20_compile.html#autotoc_md3", null ],
          [ "Option 2 (snapshot)", "page20_compile.html#autotoc_md4", null ],
          [ "Option 3 (released)", "page20_compile.html#autotoc_md5", null ],
          [ "Option 4 (binary)", "page20_compile.html#autotoc_md6", null ]
        ] ],
        [ "Compile the Code", "page20_compile.html#autotoc_md7", null ],
        [ "Compiler options (Optional)", "page20_compile.html#autotoc_md8", null ],
        [ "GPU acceleration (Optional)", "page20_compile.html#autotoc_md9", null ],
        [ "Matlab Mex Routines (Optional)", "page20_compile.html#autotoc_md10", null ],
        [ "Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon (Optional)", "page20_compile.html#autotoc_md11", null ]
      ] ],
      [ "Graphical User Interface", "page20_compile.html#autotoc_md12", [
        [ "Monitor", "page20_compile.html#autotoc_md13", null ],
        [ "Drawdaemon", "page20_compile.html#autotoc_md14", null ],
        [ "Plotting results", "page20_compile.html#autotoc_md15", null ]
      ] ],
      [ "Python Scripts", "page20_compile.html#autotoc_md16", [
        [ "Interface to MAOS", "page20_compile.html#autotoc_md17", null ]
      ] ]
    ] ],
    [ "Run simulations", "page30_run.html", [
      [ "Usage", "page30_run.html#sect-run", null ],
      [ "Configuration Files", "page30_run.html#sect-config", null ],
      [ "Sample Runs", "page30_run.html#sect-exe", null ],
      [ "Advanced configuration", "page30_run.html#advanced", [
        [ "Specifying Surface OPDs", "page30_run.html#sect-surface", null ],
        [ "WFS Configuration", "page30_run.html#sect-wfs", null ],
        [ "Point Spread Function", "page30_run.html#sect-perfevl", null ],
        [ "Actuator Slaving", "page30_run.html#sect-act", null ],
        [ "Sodium range variation", "page30_run.html#sect-sodium", null ]
      ] ],
      [ "Sky coverage", "page30_run.html#skycoverage", null ]
    ] ],
    [ "Examples", "page33_example.html", "page33_example" ],
    [ "Simulation Results", "page40_results.html", [
      [ "RMS WFE", "page40_results.html#maosres", [
        [ "Python", "page40_results.html#sect-python", null ],
        [ "IDL", "page40_results.html#sect-idl", null ],
        [ "FITS", "page40_results.html#sect-fits", null ],
        [ "Python", "page40_results.html#autotoc_md19", null ],
        [ "Matlab", "page40_results.html#autotoc_md20", null ]
      ] ],
      [ "Plotting Results", "page40_results.html#autotoc_md21", null ],
      [ "Reading Results", "page40_results.html#autotoc_md22", [
        [ ".bin file format", "page40_results.html#autotoc_md23", null ],
        [ "MATLAB", "page40_results.html#autotoc_md24", null ]
      ] ],
      [ "Result Files", "page40_results.html#autotoc_md25", [
        [ "Wavefront error", "page40_results.html#autotoc_md26", null ],
        [ "Split tomography", "page40_results.html#autotoc_md27", null ],
        [ "Log files", "page40_results.html#autotoc_md28", null ],
        [ "PSF", "page40_results.html#autotoc_md29", null ],
        [ "Other", "page40_results.html#autotoc_md30", null ]
      ] ],
      [ "Geometry Data", "page40_results.html#geometry", null ],
      [ "Telemetry Data", "page40_results.html#telemetry", null ]
    ] ],
    [ "NFIRAOS Performance", "page43_nfiraos.html", [
      [ "Turbulence profile", "page43_nfiraos.html#autotoc_md31", null ],
      [ "NFIRAOS AO Performance", "page43_nfiraos.html#autotoc_md32", [
        [ "NFIRAOS", "page43_nfiraos.html#autotoc_md33", null ],
        [ "IRIS Imager", "page43_nfiraos.html#autotoc_md34", [
          [ "Wavefront error", "page43_nfiraos.html#autotoc_md35", null ],
          [ "J band Strehl Ratio", "page43_nfiraos.html#autotoc_md36", null ],
          [ "H band Strehl Ratio", "page43_nfiraos.html#autotoc_md37", null ],
          [ "K band Strehl Ratio", "page43_nfiraos.html#autotoc_md38", null ]
        ] ],
        [ "MODHIS", "page43_nfiraos.html#autotoc_md39", null ]
      ] ],
      [ "NFIRAOS+", "page43_nfiraos.html#autotoc_md40", null ]
    ] ],
    [ "Algorithms", "algorithm.html", [
      [ "DM Actuator Influence Function", "algorithm.html#sect-dm-actuator", [
        [ "Linear influence function", "algorithm.html#autotoc_md41", null ],
        [ "Cubic influence function", "algorithm.html#autotoc_md42", null ]
      ] ],
      [ "DM Hysteresis", "algorithm.html#hysteresis", null ],
      [ "Physical Optics Beam Propagation", "algorithm.html#pop", [
        [ "Maxwell Equation", "algorithm.html#autotoc_md43", null ],
        [ "Fresnel diffraction integral", "algorithm.html#autotoc_md44", null ],
        [ "Fresnel approximation", "algorithm.html#autotoc_md45", [
          [ "Angular Spectrum", "algorithm.html#autotoc_md46", null ],
          [ "Single FFT", "algorithm.html#autotoc_md47", null ]
        ] ],
        [ "Fraunhofer approximation", "algorithm.html#autotoc_md48", null ],
        [ "Sphere to sphere propagation", "algorithm.html#autotoc_md49", null ]
      ] ]
    ] ],
    [ "Development", "page90_devel.html", "page90_devel" ],
    [ "Todo List", "todo.html", null ],
    [ "Data Structures", "annotated.html", [
      [ "Data Structures", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Data Fields", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", null ],
        [ "Variables", "functions_vars.html", "functions_vars" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "Globals", "globals.html", [
        [ "All", "globals.html", "globals_dup" ],
        [ "Functions", "globals_func.html", "globals_func" ],
        [ "Variables", "globals_vars.html", null ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"ahst_8h.html",
"classcgtmp__t.html#a02931b311da7ea7289c50cc40dc5b1b5",
"classcurecon__t.html#a9bd713d23ad6f0f485d19c81c0aa569e",
"cmath_8h.html#a4a8db9c35e106c64a679ae96beb61b66",
"cublas_8h.html#aa8395302f8567d745bfa7a5137eb6a4f",
"dmath_8h.html#ae490be1308b3596d8f88b8acbf3f4bce",
"lib_2accphi_8h.html#a7957c99547a5a45a7fb7418c68e837e1",
"maos_2types_8h.html#a0444e9675ebda314d3fc589a77caf342",
"maos_2types_8h.html#a8131af7aebef76af45f891fadc512ab1",
"maos_2types_8h.html#afd3bf7999f66779351e1ce2b60aeb42f",
"muv_8h.html#ad69fad8e7ca6f125c5697622dcb378bd",
"parms_8h.html#a2f9d6fdac75bb33eec9d79c822d522a7",
"parms_8h.html#a920d6ec25baf5f2b55cdfe1216a9363f",
"parms_8h.html#aed26032ff7e1dc04458663822ca878fe",
"recon__utils_8h.html#af7b6fd38efb1954ea36cc7b9082b126b",
"smath_8h.html#a8e7e23073d6106378e0e7ddb687f249a",
"surf_8h.html",
"type_8h.html#a21de26caa6bcfc936378c4e45d235bd9",
"type_8h.html#a9473db81d7af66f5e2a8319424834877",
"type_8h.html#structsspcell_8____unnamed_____8____unnamed____",
"zmath_8h.html#aa6255107cdd0413d115fe16632eb8cda"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';