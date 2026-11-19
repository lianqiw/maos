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
    [ "Overview", "index.html", null ],
    [ "CLAUDE.md", "d6/d23/md__c_l_a_u_d_e.html", [
      [ "Build", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md1", null ],
      [ "Tests", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md2", null ],
      [ "Running", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md3", null ],
      [ "Architecture", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md4", [
        [ "Configuration parameters", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md5", null ]
      ] ],
      [ "Conventions (from <span class=\"tt\">docs/91_basics.md</span>)", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md6", null ]
    ] ],
    [ "Possible bugs reported (code unchanged)", "d3/da0/md_claude__bugs__20261119.html", [
      [ "earlier rounds", "d3/da0/md_claude__bugs__20261119.html#autotoc_md8", null ],
      [ "batch15 tools/scheduler*", "d3/da0/md_claude__bugs__20261119.html#autotoc_md9", null ],
      [ "batch11 sys/mex/skyc", "d3/da0/md_claude__bugs__20261119.html#autotoc_md10", null ],
      [ "batch10 math/sp*, zfarr, mex/interface.h", "d3/da0/md_claude__bugs__20261119.html#autotoc_md11", null ],
      [ "batch12 sys/mem,misc,path,process*,readstr,scheduler_client,sha1", "d3/da0/md_claude__bugs__20261119.html#autotoc_md12", null ],
      [ "batch08 math/loc,map,mat,matbin,matcomp", "d3/da0/md_claude__bugs__20261119.html#autotoc_md13", null ],
      [ "batch14 tools (monitor, ethtest, memdbg, mvm_cpu, psfgc...)", "d3/da0/md_claude__bugs__20261119.html#autotoc_md14", null ],
      [ "batch13 sys/sock*,thread*,zstdio, tools/drawdaemon*", "d3/da0/md_claude__bugs__20261119.html#autotoc_md15", null ],
      [ "batch09 math/mathmisc,numtype,random,matmmap,matmath", "d3/da0/md_claude__bugs__20261119.html#autotoc_md16", null ],
      [ "batch07 maos/sim*,surf,twfs,utils,wfsgrad,wfsints, math/amoeba,bessik,blas,cell,fft", "d3/da0/md_claude__bugs__20261119.html#autotoc_md17", null ],
      [ "batch05 lib/servo,slaving,stfun,turbulence,zernike; maos/ahst,dither_utils,filter,fdpcg,maos,moao,mvm_client", "d3/da0/md_claude__bugs__20261119.html#autotoc_md18", null ],
      [ "batch02 cuda/recon, cuda/sim, cuda/test (part)", "d3/da0/md_claude__bugs__20261119.html#autotoc_md19", null ],
      [ "batch03 cuda/test, lib/accphi..libmisc", "d3/da0/md_claude__bugs__20261119.html#autotoc_md20", null ],
      [ "batch01 cuda/math", "d3/da0/md_claude__bugs__20261119.html#autotoc_md21", null ],
      [ "batch04 lib/locfft..readcfg", "d3/da0/md_claude__bugs__20261119.html#autotoc_md22", null ],
      [ "batch06 maos/parms..setup", "d3/da0/md_claude__bugs__20261119.html#autotoc_md23", null ]
    ] ],
    [ "MAOS bug-fix proposals (for review)", "df/d8d/md_claude__review__20261119.html", [
      [ "Files", "df/d8d/md_claude__review__20261119.html#autotoc_md25", null ],
      [ "Verification of the combined patch", "df/d8d/md_claude__review__20261119.html#autotoc_md26", null ],
      [ "Counts", "df/d8d/md_claude__review__20261119.html#autotoc_md27", null ],
      [ "Your earlier fixes in 88beaf2da that were incomplete", "df/d8d/md_claude__review__20261119.html#autotoc_md28", null ],
      [ "Fixes that can change simulation results", "df/d8d/md_claude__review__20261119.html#autotoc_md29", null ],
      [ "Security", "df/d8d/md_claude__review__20261119.html#autotoc_md30", null ],
      [ "Needs your decision", "df/d8d/md_claude__review__20261119.html#autotoc_md31", null ],
      [ "Noticed but not changed", "df/d8d/md_claude__review__20261119.html#autotoc_md32", null ],
      [ "Per-group verdict tables", "df/d8d/md_claude__review__20261119.html#autotoc_md34", [
        [ "Group: sys", "df/d8d/md_claude__review__20261119.html#autotoc_md35", [
          [ "sys group: fixes", "df/d8d/md_claude__review__20261119.html#autotoc_md36", null ]
        ] ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md37", null ],
        [ "Group: math", "df/d8d/md_claude__review__20261119.html#autotoc_md38", [
          [ "math/ fixes", "df/d8d/md_claude__review__20261119.html#autotoc_md39", null ]
        ] ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md40", null ],
        [ "Related notes (outside math/, not changed)", "df/d8d/md_claude__review__20261119.html#autotoc_md41", null ],
        [ "Group: lib", "df/d8d/md_claude__review__20261119.html#autotoc_md42", [
          [ "lib group: bug triage", "df/d8d/md_claude__review__20261119.html#autotoc_md43", null ]
        ] ],
        [ "mkh.c renormalization check (requested)", "df/d8d/md_claude__review__20261119.html#autotoc_md44", null ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md45", null ],
        [ "Group: maos_skyc", "df/d8d/md_claude__review__20261119.html#autotoc_md46", [
          [ "maos_skyc fixes", "df/d8d/md_claude__review__20261119.html#autotoc_md47", null ]
        ] ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md48", null ],
        [ "Group: cuda", "df/d8d/md_claude__review__20261119.html#autotoc_md49", [
          [ "cuda group: bug review", "df/d8d/md_claude__review__20261119.html#autotoc_md50", null ]
        ] ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md51", null ],
        [ "High-impact items", "df/d8d/md_claude__review__20261119.html#autotoc_md52", null ],
        [ "Group: tools_mex", "df/d8d/md_claude__review__20261119.html#autotoc_md53", [
          [ "tools_mex fixes (base: 88beaf2da)", "df/d8d/md_claude__review__20261119.html#autotoc_md54", null ]
        ] ],
        [ "Needs decision", "df/d8d/md_claude__review__20261119.html#autotoc_md55", null ],
        [ "Build check", "df/d8d/md_claude__review__20261119.html#autotoc_md56", null ]
      ] ]
    ] ],
    [ "Introduction", "da/d9f/page10_intro.html", null ],
    [ "Compile the Code", "d3/d4c/page20_compile.html", [
      [ "Compiling", "d3/d4c/page20_compile.html#autotoc_md57", [
        [ "Prerequisites", "d3/d4c/page20_compile.html#autotoc_md58", null ],
        [ "Download the code", "d3/d4c/page20_compile.html#autotoc_md59", [
          [ "Option 1 (preferred)", "d3/d4c/page20_compile.html#autotoc_md60", null ],
          [ "Option 2 (snapshot)", "d3/d4c/page20_compile.html#autotoc_md61", null ],
          [ "Option 3 (released)", "d3/d4c/page20_compile.html#autotoc_md62", null ],
          [ "Option 4 (binary)", "d3/d4c/page20_compile.html#autotoc_md63", null ]
        ] ],
        [ "Compile the Code", "d3/d4c/page20_compile.html#autotoc_md64", null ],
        [ "Compiler options (Optional)", "d3/d4c/page20_compile.html#autotoc_md65", null ],
        [ "GPU acceleration (Optional)", "d3/d4c/page20_compile.html#autotoc_md66", null ],
        [ "Matlab Mex Routines (Optional)", "d3/d4c/page20_compile.html#autotoc_md67", null ],
        [ "Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon (Optional)", "d3/d4c/page20_compile.html#autotoc_md68", null ]
      ] ],
      [ "Browser User Interface", "d3/d4c/page20_compile.html#autotoc_md69", null ],
      [ "Graphical User Interface", "d3/d4c/page20_compile.html#autotoc_md70", [
        [ "Monitor", "d3/d4c/page20_compile.html#autotoc_md71", null ],
        [ "Drawdaemon", "d3/d4c/page20_compile.html#autotoc_md72", null ],
        [ "Plotting results", "d3/d4c/page20_compile.html#autotoc_md73", null ]
      ] ],
      [ "Python Scripts", "d3/d4c/page20_compile.html#autotoc_md74", [
        [ "Interface to MAOS", "d3/d4c/page20_compile.html#autotoc_md75", null ]
      ] ]
    ] ],
    [ "Run simulations", "da/dc3/page30_run.html", [
      [ "Usage", "da/dc3/page30_run.html#sect-run", null ],
      [ "Configuration Files", "da/dc3/page30_run.html#sect-config", null ],
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
      [ "RMS WFE", "dd/da9/page40_results.html#autotoc_md77", [
        [ "Python", "dd/da9/page40_results.html#autotoc_md78", null ],
        [ "Matlab", "dd/da9/page40_results.html#autotoc_md79", null ]
      ] ],
      [ "Plotting Results", "dd/da9/page40_results.html#autotoc_md80", null ],
      [ "Reading Results", "dd/da9/page40_results.html#autotoc_md81", [
        [ ".bin file format", "dd/da9/page40_results.html#autotoc_md82", null ],
        [ "MATLAB", "dd/da9/page40_results.html#autotoc_md83", null ],
        [ "Python", "dd/da9/page40_results.html#sect-python", null ],
        [ "IDL", "dd/da9/page40_results.html#sect-idl", null ],
        [ "FITS", "dd/da9/page40_results.html#sect-fits", null ]
      ] ],
      [ "Result Files", "dd/da9/page40_results.html#autotoc_md84", [
        [ "Wavefront error", "dd/da9/page40_results.html#autotoc_md85", null ],
        [ "Split tomography results", "dd/da9/page40_results.html#autotoc_md86", null ],
        [ "Log files", "dd/da9/page40_results.html#autotoc_md87", null ],
        [ "PSF", "dd/da9/page40_results.html#autotoc_md88", null ],
        [ "Subaperture Pixels", "dd/da9/page40_results.html#autotoc_md89", null ],
        [ "Other", "dd/da9/page40_results.html#autotoc_md90", null ]
      ] ],
      [ "Geometry Data", "dd/da9/page40_results.html#geometry", null ],
      [ "Telemetry Data", "dd/da9/page40_results.html#telemetry", null ]
    ] ],
    [ "NFIRAOS Performance", "d6/d3f/page43_nfiraos.html", [
      [ "Turbulence profile", "d6/d3f/page43_nfiraos.html#autotoc_md91", null ],
      [ "NFIRAOS AO Performance", "d6/d3f/page43_nfiraos.html#autotoc_md92", [
        [ "NFIRAOS", "d6/d3f/page43_nfiraos.html#autotoc_md93", null ],
        [ "IRIS Imager", "d6/d3f/page43_nfiraos.html#autotoc_md94", [
          [ "Wavefront error", "d6/d3f/page43_nfiraos.html#autotoc_md95", null ],
          [ "J band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md96", null ],
          [ "H band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md97", null ],
          [ "K band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md98", null ]
        ] ],
        [ "MODHIS", "d6/d3f/page43_nfiraos.html#autotoc_md99", null ]
      ] ],
      [ "NFIRAOS+", "d6/d3f/page43_nfiraos.html#autotoc_md100", null ]
    ] ],
    [ "Algorithms", "d4/d15/algorithm.html", [
      [ "DM Actuator Influence Function", "d4/d15/algorithm.html#sect-dm-actuator", [
        [ "Linear influence function", "d4/d15/algorithm.html#autotoc_md101", null ],
        [ "Cubic influence function", "d4/d15/algorithm.html#autotoc_md102", null ]
      ] ],
      [ "DM Hysteresis", "d4/d15/algorithm.html#hysteresis", null ],
      [ "Physical Optics Beam Propagation", "d4/d15/algorithm.html#pop", [
        [ "Maxwell Equation", "d4/d15/algorithm.html#autotoc_md103", null ],
        [ "Fresnel diffraction integral", "d4/d15/algorithm.html#autotoc_md104", null ],
        [ "Fresnel approximation", "d4/d15/algorithm.html#autotoc_md105", [
          [ "Angular Spectrum", "d4/d15/algorithm.html#autotoc_md106", null ],
          [ "Single FFT", "d4/d15/algorithm.html#autotoc_md107", null ]
        ] ],
        [ "Fraunhofer approximation", "d4/d15/algorithm.html#autotoc_md108", null ],
        [ "Sphere to sphere propagation", "d4/d15/algorithm.html#autotoc_md109", null ]
      ] ]
    ] ],
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
"d2/d05/powfs_8h.html#a1a9c34e6fc651433c73473e0e039dc63",
"d3/dc7/muv_8h.html#a680b2b553061f7be6f8925b41ed9260c",
"d4/d38/parms_8h.html#a1618076f04cfcee3d65c29f0451b0252",
"d4/d38/parms_8h.html#a7e5fb4451f9c4ee844e8013b0aefc0a0",
"d4/d38/parms_8h.html#ae093c8027457ac84653fbf51f5e325e6",
"d5/d98/page91_basics.html#sect-guide",
"d6/d94/maos_2types_8h.html#a5d0b1fa0a3a88bb8532ccd7631b3291f",
"d6/d94/maos_2types_8h.html#ae69c496d45aeef46c543e0306592f029",
"da/d0f/keck_base.html",
"da/dfb/dmath_8h.html#ad6518403ee2ca0e2d52eaf2a0813e6c7",
"dc/db7/map_8h.html",
"de/d07/cmath_8h.html#a662b273b76b56a002d9bdf6b5c6086ee",
"de/de7/smath_8h.html#a1cecf05e2bbf08c2ded9d41c27204cfe",
"df/d2b/mathmisc_8h.html#a1d7e5986f277a0e9546311fef2526414",
"df/def/type_8h.html#a85058f9c06c6ac49068e2b77ffd67265",
"df/def/type_8h.html#df/d5e/structlcccell"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';