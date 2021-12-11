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
      [ "Requirements", "page20_compile.html#autotoc_md0", null ],
      [ "Preparing the folders and compiling", "page20_compile.html#autotoc_md1", null ],
      [ "Graphical User Interface", "page20_compile.html#autotoc_md2", [
        [ "Compiler options", "page20_compile.html#sect-compiler", null ],
        [ "GPU acceleration", "page20_compile.html#sect-cuda", null ],
        [ "Matlab Mex Routines", "page20_compile.html#sect-mex", null ],
        [ "Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon", "page20_compile.html#sect-mac-gtk", null ],
        [ "Monitor", "page20_compile.html#autotoc_md3", null ],
        [ "Drawdaemon", "page20_compile.html#autotoc_md4", null ],
        [ "Plotting results", "page20_compile.html#autotoc_md5", null ]
      ] ],
      [ "Python Scripts", "page20_compile.html#autotoc_md6", [
        [ "Interface to MAOS", "page20_compile.html#autotoc_md7", null ]
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
      ] ]
    ] ],
    [ "Examples", "page33_example.html", "page33_example" ],
    [ "Simulation Results", "page40_results.html", [
      [ "RMS WFE", "page40_results.html#maosres", [
        [ "Python", "page40_results.html#sect-python", null ],
        [ "IDL", "page40_results.html#sect-idl", null ],
        [ "FITS", "page40_results.html#sect-fits", null ],
        [ "Python", "page40_results.html#autotoc_md9", null ],
        [ "Matlab", "page40_results.html#autotoc_md10", null ]
      ] ],
      [ "Plotting Results", "page40_results.html#autotoc_md11", null ],
      [ "Reading Results", "page40_results.html#autotoc_md12", [
        [ ".bin file format", "page40_results.html#autotoc_md13", null ],
        [ "MATLAB", "page40_results.html#autotoc_md14", null ]
      ] ],
      [ "Result Files", "page40_results.html#autotoc_md15", [
        [ "Wavefront error", "page40_results.html#autotoc_md16", null ],
        [ "Split tomography", "page40_results.html#autotoc_md17", null ],
        [ "Log files", "page40_results.html#autotoc_md18", null ],
        [ "PSF", "page40_results.html#autotoc_md19", null ],
        [ "Other", "page40_results.html#autotoc_md20", null ]
      ] ],
      [ "Geometry Data", "page40_results.html#geometry", null ],
      [ "Telemetry Data", "page40_results.html#telemetry", null ]
    ] ],
    [ "NFIRAOS Performance", "page43_nfiraos.html", [
      [ "Turbulence profile", "page43_nfiraos.html#autotoc_md21", null ],
      [ "AO Performance", "page43_nfiraos.html#autotoc_md22", [
        [ "NFIRAOS", "page43_nfiraos.html#autotoc_md23", null ],
        [ "IRIS Imager", "page43_nfiraos.html#autotoc_md24", [
          [ "Wavefront error", "page43_nfiraos.html#autotoc_md25", null ],
          [ "J band Strehl Ratio", "page43_nfiraos.html#autotoc_md26", null ],
          [ "H band Strehl Ratio", "page43_nfiraos.html#autotoc_md27", null ],
          [ "K band Strehl Ratio", "page43_nfiraos.html#autotoc_md28", null ]
        ] ],
        [ "MODHIS", "page43_nfiraos.html#autotoc_md29", null ]
      ] ]
    ] ],
    [ "Algorithms", "algorithm.html", [
      [ "DM Actuator Influence Function", "algorithm.html#sect-dm-actuator", [
        [ "Linear influence function", "algorithm.html#autotoc_md30", null ],
        [ "Cubic influence function", "algorithm.html#autotoc_md31", null ]
      ] ],
      [ "DM Hysteresis", "algorithm.html#hysteresis", null ],
      [ "Physical Optics Beam Propagation", "algorithm.html#pop", [
        [ "Maxwell Equation", "algorithm.html#autotoc_md32", null ],
        [ "Fresnel diffraction integral", "algorithm.html#autotoc_md33", null ],
        [ "Fresnel approximation", "algorithm.html#autotoc_md34", [
          [ "Angular Spectrum", "algorithm.html#autotoc_md35", null ],
          [ "Single FFT", "algorithm.html#autotoc_md36", null ]
        ] ],
        [ "Fraunhofer approximation", "algorithm.html#autotoc_md37", null ],
        [ "Sphere to sphere propagation", "algorithm.html#autotoc_md38", null ]
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
".html",
"classcuda__recon_1_1cufit__grid.html#aa2898dcc7210a5e3a1cd708b5c7d8acf",
"classculoc__t.html",
"cmath_8h.html#a69f23aabe8d7a67eb05bb6699512517a",
"dmath_8h.html#a10830fb41f473ae929a5aee34361dec1",
"draw_8h.html#af6a614e128ea550040fe980a1b1922e6",
"lmath_8h.html#a15ceffb00e7bb7e57793410cecd33cbd",
"maos_2types_8h.html#a28d8f689a09cdfaa209bae03ee8d987d",
"maos_2types_8h.html#aaa286b3af6f54cc7c0157956a85e90c3",
"mathdef_8h.html#a305fa5f76f2bb165561515f9971475a9",
"mvm__trans_8cu.html#a89e6c5c637455fc0b9c1b43cb1d0e78e",
"page43_nfiraos.html#autotoc_md25",
"parms_8h.html#a6172520813e1bd7d0a7464cb0c6dd3fa",
"parms_8h.html#ac08de94a88488ca5f106a8db93990bc4",
"pywfs_8h.html#a268491858839c93d7e3f79ed514f3b47",
"slaving_8h.html#a7bbac27715f60f0b22b05b548ed840e5",
"smath_8h.html#af8f3f1729491c128bbccf8e3b9eba724",
"test_8cu.html",
"type_8h.html#a546c88ba8a0e1904479cc349aef10feb",
"type_8h.html#ab80bb7740288fda1f201890375a60c8f",
"zmath_8h.html#a29d71d5b8f782fca890c372a12e4e50b"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';