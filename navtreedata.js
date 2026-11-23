/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "MAOS", "index.html", [
    [ "Overview", "index.html", "index" ],
    [ "Introduction", "da/d9f/page10_intro.html", [
      [ "What is modeled", "da/d9f/page10_intro.html#intro-model", [
        [ "Configuration", "da/d9f/page10_intro.html#intro-config", null ],
        [ "Atmosphere", "da/d9f/page10_intro.html#intro-atm", null ],
        [ "Telescope", "da/d9f/page10_intro.html#intro-telescope", null ],
        [ "Wavefront sensors", "da/d9f/page10_intro.html#intro-wfs", null ],
        [ "Deformable mirrors", "da/d9f/page10_intro.html#intro-dm", null ],
        [ "Reconstruction and control", "da/d9f/page10_intro.html#intro-recon", null ],
        [ "Performance evaluation", "da/d9f/page10_intro.html#intro-perf", null ]
      ] ],
      [ "Implementation", "da/d9f/page10_intro.html#intro-impl", null ],
      [ "Source tree", "da/d9f/page10_intro.html#intro-source", null ]
    ] ],
    [ "Compile the Code", "d3/d4c/page20_compile.html", [
      [ "Quick start", "d3/d4c/page20_compile.html#compile-quick", null ],
      [ "Prerequisites", "d3/d4c/page20_compile.html#compile-prereq", null ],
      [ "Download the code", "d3/d4c/page20_compile.html#compile-download", [
        [ "Git (preferred)", "d3/d4c/page20_compile.html#compile-git", null ],
        [ "Tagged snapshot", "d3/d4c/page20_compile.html#compile-snapshot", null ],
        [ "Release tarball", "d3/d4c/page20_compile.html#compile-release", null ],
        [ "Pre-compiled binaries", "d3/d4c/page20_compile.html#compile-binary", null ]
      ] ],
      [ "Build", "d3/d4c/page20_compile.html#compile-build", null ],
      [ "Configure options", "d3/d4c/page20_compile.html#compile-options", [
        [ "Compilers", "d3/d4c/page20_compile.html#compile-compilers", null ],
        [ "C++ build", "d3/d4c/page20_compile.html#compile-cpp", null ],
        [ "Single precision", "d3/d4c/page20_compile.html#compile-single", null ]
      ] ],
      [ "External libraries", "d3/d4c/page20_compile.html#compile-deps", null ],
      [ "GPU acceleration", "d3/d4c/page20_compile.html#compile-gpu", null ],
      [ "MATLAB mex interface", "d3/d4c/page20_compile.html#compile-mex", null ],
      [ "Unit tests and documentation", "d3/d4c/page20_compile.html#compile-test", null ],
      [ "Continuous integration", "d3/d4c/page20_compile.html#compile-ci", null ],
      [ "Troubleshooting", "d3/d4c/page20_compile.html#compile-trouble", null ],
      [ "Browser based monitor", "d3/d4c/page20_compile.html#compile-browser", null ],
      [ "Graphical user interface", "d3/d4c/page20_compile.html#compile-gui", [
        [ "Monitor", "d3/d4c/page20_compile.html#compile-monitor", null ],
        [ "Drawdaemon", "d3/d4c/page20_compile.html#compile-drawdaemon", null ],
        [ "Plotting results", "d3/d4c/page20_compile.html#compile-drawres", null ]
      ] ],
      [ "Python scripts", "d3/d4c/page20_compile.html#compile-python", [
        [ "Interface to a running maos", "d3/d4c/page20_compile.html#compile-maos-client", null ]
      ] ]
    ] ],
    [ "Run simulations", "da/dc3/page30_run.html", [
      [ "Usage", "da/dc3/page30_run.html#sect-run", [
        [ "Command line options", "da/dc3/page30_run.html#sect-run-options", null ],
        [ "Scheduler and detached runs", "da/dc3/page30_run.html#sect-run-sched", null ],
        [ "GPU selection", "da/dc3/page30_run.html#sect-run-gpu", null ],
        [ "Seeds, restarting, and multiple runs", "da/dc3/page30_run.html#sect-run-seeds", null ],
        [ "Output folder", "da/dc3/page30_run.html#sect-run-output", null ],
        [ "Environment variables", "da/dc3/page30_run.html#sect-run-env", null ]
      ] ],
      [ "Configuration Files", "da/dc3/page30_run.html#sect-config", [
        [ "Syntax", "da/dc3/page30_run.html#sect-config-syntax", null ],
        [ "Overriding and checking keys", "da/dc3/page30_run.html#sect-config-override", null ],
        [ "Search path", "da/dc3/page30_run.html#sect-config-path", null ],
        [ "Organization", "da/dc3/page30_run.html#sect-config-org", null ],
        [ "Guide to the main parameters", "da/dc3/page30_run.html#sect-config-params", null ]
      ] ],
      [ "Sample Runs", "da/dc3/page30_run.html#sect-exe", null ],
      [ "Error Budget Breakdown", "da/dc3/page30_run.html#sect-ebb", null ],
      [ "Advanced configuration", "da/dc3/page30_run.html#advanced", [
        [ "Specifying Surface OPDs", "da/dc3/page30_run.html#sect-surface", null ],
        [ "WFS Configuration", "da/dc3/page30_run.html#sect-wfs", null ],
        [ "Point Spread Function", "da/dc3/page30_run.html#sect-perfevl", null ],
        [ "Actuator Slaving", "da/dc3/page30_run.html#sect-act", null ],
        [ "Sodium range variation", "da/dc3/page30_run.html#sect-sodium", null ],
        [ "Rayleigh backscatter", "da/dc3/page30_run.html#pixel-background", null ]
      ] ],
      [ "Sky coverage", "da/dc3/page30_run.html#maos_skycoverage", null ]
    ] ],
    [ "Examples", "d5/d29/page33_example.html", "d5/d29/page33_example" ],
    [ "Simulation Results", "dd/da9/page40_results.html", [
      [ "RMS WFE", "dd/da9/page40_results.html#autotoc_md0", [
        [ "Python", "dd/da9/page40_results.html#autotoc_md1", null ],
        [ "MATLAB", "dd/da9/page40_results.html#autotoc_md2", null ]
      ] ],
      [ "Plotting Results", "dd/da9/page40_results.html#autotoc_md3", null ],
      [ "Reading Results", "dd/da9/page40_results.html#autotoc_md4", [
        [ ".bin file format", "dd/da9/page40_results.html#autotoc_md5", null ],
        [ "MATLAB", "dd/da9/page40_results.html#sect-mex", null ],
        [ "Python", "dd/da9/page40_results.html#sect-python", null ],
        [ "IDL", "dd/da9/page40_results.html#sect-idl", null ],
        [ "FITS", "dd/da9/page40_results.html#sect-fits", null ]
      ] ],
      [ "Result Files", "dd/da9/page40_results.html#autotoc_md6", [
        [ "Wavefront error", "dd/da9/page40_results.html#sect-res", null ],
        [ "Split tomography results", "dd/da9/page40_results.html#autotoc_md7", null ],
        [ "Log files", "dd/da9/page40_results.html#autotoc_md8", null ],
        [ "PSF", "dd/da9/page40_results.html#autotoc_md9", null ],
        [ "Subaperture Pixels", "dd/da9/page40_results.html#autotoc_md10", null ],
        [ "Other", "dd/da9/page40_results.html#autotoc_md11", null ]
      ] ],
      [ "Geometry Data", "dd/da9/page40_results.html#geometry", null ],
      [ "Telemetry Data", "dd/da9/page40_results.html#telemetry", [
        [ "Summary of save options", "dd/da9/page40_results.html#autotoc_md12", null ]
      ] ]
    ] ],
    [ "NFIRAOS Performance", "d6/d3f/page43_nfiraos.html", [
      [ "Turbulence profile", "d6/d3f/page43_nfiraos.html#autotoc_md13", null ],
      [ "NFIRAOS AO Performance", "d6/d3f/page43_nfiraos.html#autotoc_md14", [
        [ "NFIRAOS", "d6/d3f/page43_nfiraos.html#autotoc_md15", null ],
        [ "IRIS Imager", "d6/d3f/page43_nfiraos.html#autotoc_md16", [
          [ "Wavefront error", "d6/d3f/page43_nfiraos.html#autotoc_md17", null ],
          [ "J band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md18", null ],
          [ "H band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md19", null ],
          [ "K band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md20", null ]
        ] ],
        [ "MODHIS", "d6/d3f/page43_nfiraos.html#autotoc_md21", null ]
      ] ],
      [ "NFIRAOS+", "d6/d3f/page43_nfiraos.html#autotoc_md22", null ]
    ] ],
    [ "Algorithms", "d4/d15/algorithm.html", "d4/d15/algorithm" ],
    [ "Sky Coverage", "d5/d4c/skycoverage.html", [
      [ "MAOS presimulation", "d5/d4c/skycoverage.html#presim", null ],
      [ "Sky coverage post processing simulation", "d5/d4c/skycoverage.html#postsim", null ],
      [ "Sky coverage results", "d5/d4c/skycoverage.html#skyres", null ]
    ] ],
    [ "Development", "de/d19/page90_devel.html", "de/d19/page90_devel" ],
    [ "README", "d3/dcc/md__r_e_a_d_m_e.html", null ],
    [ "Todo List", "dd/da0/todo.html", null ],
    [ "Data Structures", "annotated.html", [
      [ "Data Structures", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Data Fields", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Related Symbols", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "Globals", "globals.html", [
        [ "All", "globals.html", "globals_dup" ],
        [ "Functions", "globals_func.html", "globals_func" ],
        [ "Variables", "globals_vars.html", null ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Enumerations", "globals_enum.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"d1/d83/loc_8h.html#ae3ff7dce2c0f49a6e1afb12d559fa4a3",
"d3/d36/zmath_8h.html#a4ec8d6287c393bf17360a5f53f04c762",
"d3/dbe/evl_x.html",
"d4/d38/parms_8h.html#a1310f644d109168bdc986832854290e2",
"d4/d38/parms_8h.html#a7bf4ca9968bca8ce8b4bbb3aa41d7565",
"d4/d38/parms_8h.html#adf33bd76b9ccf6fbb5cdaa7d6ae6641a",
"d5/d85/surf_8h.html#ab9d628f58bc2c7aa07274f9f00462107",
"d6/d94/maos_2types_8h.html#a5d0b1fa0a3a88bb8532ccd7631b3291f",
"d6/d94/maos_2types_8h.html#aeb5b9028e45fb2109c6d262675a1ea93",
"d9/d6e/page96_unittest.html",
"da/dfb/dmath_8h.html#a7ff47722a9401a5bf0e98b3be7493920",
"dc/d59/cn2est_8h.html#a742aa5a993650df584a3c505336266d5",
"dd/dc1/lmath_8h.html#a5b9920077ade8f512d7c91820a5c1872",
"de/d07/cmath_8h.html#aeb46f604c51f21d2ae0de2f764f3f14a",
"de/de7/smath_8h.html#a9fec2042f28f1838095d50941f1e7e43",
"df/def/type_8h.html#a42e2503f68240fc7269809cf7e7400cf",
"df/def/type_8h.html#d6/de5/structcccell"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';