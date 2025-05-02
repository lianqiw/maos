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
    [ "Sky Coverage", "skycoverage.html", [
      [ "presimulation", "skycoverage.html#MAOS", null ],
      [ "coverage post processing simulation", "skycoverage.html#Sky", null ]
    ] ],
    [ "Development", "page90_devel.html", "page90_devel" ],
    [ "README", "md__r_e_a_d_m_e.html", null ],
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
"class_num_cell.html#a611a5f54fa71145ce20db1ea85fa9156",
"classcurecon__geom.html",
"cmath_8h.html#a1f6ecbeecb56f05d855f3e1f21342dc2",
"cn2est_8h.html#a48f8636a584276be0b5dad3b2b08dbf5",
"dmath_8h.html#abf95f69d34879322190d3ea63c73ca7d",
"kalman_8h.html#a4d3a82244fe94108b9ab5b1c57d06f38",
"locfft_8h.html#a3ed16688589c043dbe270cf928071113",
"maos_2types_8h.html#a57579226ab2934985993550e03c153c2",
"maos_2types_8h.html#ad7c3387543ce591091d91d3f894158ab",
"mcao_ngs.html",
"page40_results.html#geometry",
"parms_8h.html#a5a0944fd4c266098d168c34db1d1e8c1aaaeb75f8ad3790c0d26cc2a0dcf388f4",
"parms_8h.html#ab4d4972861877b925caec10fd38f528b",
"perfevl_8c.html#a799c7105d3e2ffd9cbcd7de50fbd16e7",
"scheduler__client_8h.html#ad67c7f31687c8712e4f6359f06c6267f",
"smath_8h.html#aa9c3e9f38857fcbd803f59a9b0430c40",
"structlmat.html",
"type_8h.html#a0bcbd3b2fc2dd5980b3bda19902edfbf",
"type_8h.html#aa12d0a13a21ea965033a305302ed8368",
"zmath_8h.html#a07fae36f0643c5069d42214431a0569d",
"zmath_8h.html#af9c7e8323c544620e68e2ee9a9a02ba3"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';