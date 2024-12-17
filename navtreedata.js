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
      [ "Error Budget Breakdown", "page30_run.html#sect-ebb", null ],
      [ "Advanced configuration", "page30_run.html#advanced", [
        [ "Specifying Surface OPDs", "page30_run.html#sect-surface", null ],
        [ "WFS Configuration", "page30_run.html#sect-wfs", null ],
        [ "Point Spread Function", "page30_run.html#sect-perfevl", null ],
        [ "Actuator Slaving", "page30_run.html#sect-act", null ],
        [ "Sodium range variation", "page30_run.html#sect-sodium", null ],
        [ "Rayleigh backscatter", "page30_run.html#pixel-background", null ]
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
"class_num_array.html#aadcad54fb14e9a4e3b9b23580938f24c",
"classcurecon__geom.html#a35471d2823ee5c33c2e26ac07ee9757e",
"cmath_8h.html#a2281038438f8486887be78c4b3a95f60",
"cn2est_8h.html#a7b84c08f588748d9c5f4582ff809b5b2",
"dmath_8h.html#acb291ac90c6280592569f8a1e1bd20dd",
"kalman_8h.html#af46e2a92d29de6efda5353657f3da3a3",
"maos_2recon_8h.html#a4bfc73b39b67ad8fd1ff3329552e80ee",
"maos_2types_8h.html#a70068bd1ecc6478dc836162e2bfbf777",
"maos_2types_8h.html#aecfdde04f1823426a4a6bde295ae6e95",
"misc_8h.html#aa8327247b00b74fe0f22822f797f5ec4",
"parms_8h.html#a036091e39c7f9549b1116b6cfc8cc90a",
"parms_8h.html#a69588f63ff48ddd8e9b0bba6c9c86ab0",
"parms_8h.html#ac6585987a62ac18077518e4edc4559a5",
"process_8h.html#ab669b7d706b102f1f5cc460624f626b3",
"sim_8h.html#a2f3ab210392cf3ba70bcf8df2ccbfde3",
"smath_8h.html#aeee30e5c2c5f76cb95bb8f0eb895d302",
"thread_8h.html#a662929086396591a2ce1ad15f451c4e7",
"type_8h.html#a6cd3d45c39f5e6bba94a8a535172d9c1",
"type_8h.html#aeb40fac426ae3bc474f3edb78bef3844",
"zmath_8h.html#a93ffc9078e81536bb2ccbd3fc74c467a"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';