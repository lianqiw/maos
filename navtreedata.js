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
    [ "CLAUDE.md", "d6/d23/md__c_l_a_u_d_e.html", [
      [ "Build", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md1", null ],
      [ "Tests", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md2", null ],
      [ "Running", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md3", null ],
      [ "Architecture", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md4", [
        [ "Configuration parameters", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md5", null ]
      ] ],
      [ "Conventions (from <span class=\"tt\">docs/91_basics.md</span>)", "d6/d23/md__c_l_a_u_d_e.html#autotoc_md6", null ]
    ] ],
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
      [ "RMS WFE", "dd/da9/page40_results.html#autotoc_md7", [
        [ "Python", "dd/da9/page40_results.html#autotoc_md8", null ],
        [ "MATLAB", "dd/da9/page40_results.html#autotoc_md9", null ]
      ] ],
      [ "Plotting Results", "dd/da9/page40_results.html#autotoc_md10", null ],
      [ "Reading Results", "dd/da9/page40_results.html#autotoc_md11", [
        [ ".bin file format", "dd/da9/page40_results.html#autotoc_md12", null ],
        [ "MATLAB", "dd/da9/page40_results.html#sect-mex", null ],
        [ "Python", "dd/da9/page40_results.html#sect-python", null ],
        [ "IDL", "dd/da9/page40_results.html#sect-idl", null ],
        [ "FITS", "dd/da9/page40_results.html#sect-fits", null ]
      ] ],
      [ "Result Files", "dd/da9/page40_results.html#autotoc_md13", [
        [ "Wavefront error", "dd/da9/page40_results.html#sect-res", null ],
        [ "Split tomography results", "dd/da9/page40_results.html#autotoc_md14", null ],
        [ "Log files", "dd/da9/page40_results.html#autotoc_md15", null ],
        [ "PSF", "dd/da9/page40_results.html#autotoc_md16", null ],
        [ "Subaperture Pixels", "dd/da9/page40_results.html#autotoc_md17", null ],
        [ "Other", "dd/da9/page40_results.html#autotoc_md18", null ]
      ] ],
      [ "Geometry Data", "dd/da9/page40_results.html#geometry", null ],
      [ "Telemetry Data", "dd/da9/page40_results.html#telemetry", [
        [ "Summary of save options", "dd/da9/page40_results.html#autotoc_md19", null ]
      ] ]
    ] ],
    [ "NFIRAOS Performance", "d6/d3f/page43_nfiraos.html", [
      [ "Turbulence profile", "d6/d3f/page43_nfiraos.html#autotoc_md20", null ],
      [ "NFIRAOS AO Performance", "d6/d3f/page43_nfiraos.html#autotoc_md21", [
        [ "NFIRAOS", "d6/d3f/page43_nfiraos.html#autotoc_md22", null ],
        [ "IRIS Imager", "d6/d3f/page43_nfiraos.html#autotoc_md23", [
          [ "Wavefront error", "d6/d3f/page43_nfiraos.html#autotoc_md24", null ],
          [ "J band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md25", null ],
          [ "H band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md26", null ],
          [ "K band Strehl Ratio", "d6/d3f/page43_nfiraos.html#autotoc_md27", null ]
        ] ],
        [ "MODHIS", "d6/d3f/page43_nfiraos.html#autotoc_md28", null ]
      ] ],
      [ "NFIRAOS+", "d6/d3f/page43_nfiraos.html#autotoc_md29", null ]
    ] ],
    [ "Algorithms", "d4/d15/algorithm.html", "d4/d15/algorithm" ],
    [ "Sky Coverage", "d5/d4c/skycoverage.html", [
      [ "MAOS presimulation", "d5/d4c/skycoverage.html#presim", null ],
      [ "Sky coverage post processing simulation", "d5/d4c/skycoverage.html#postsim", null ],
      [ "Sky coverage results", "d5/d4c/skycoverage.html#skyres", null ]
    ] ],
    [ "Development", "de/d19/page90_devel.html", "de/d19/page90_devel" ],
    [ "Bug fixes, unit tests and documentation by Claude, 2026-11-19 to 2026-11-20", "db/d87/md_docs_2claude__changes__20261120.html", [
      [ "Verification", "db/d87/md_docs_2claude__changes__20261120.html#claude_verification", null ],
      [ "Changes that can alter results", "db/d87/md_docs_2claude__changes__20261120.html#claude_results", null ],
      [ "Not changed", "db/d87/md_docs_2claude__changes__20261120.html#claude_not_changed", null ],
      [ "First review of the reported bugs", "db/d87/md_docs_2claude__changes__20261120.html#claude_first_review", [
        [ "Verification at the time", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md34", null ],
        [ "Counts", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md35", null ],
        [ "Your earlier fixes in 88beaf2da that were incomplete", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md36", null ],
        [ "Fixes that can change simulation results", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md37", null ],
        [ "Security", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md38", null ],
        [ "Decisions (option applied unless noted)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md39", null ],
        [ "Noticed but not changed", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md40", null ],
        [ "Per-group verdict tables", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md41", null ],
        [ "Group: sys", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md42", [
          [ "sys group: fixes", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md43", null ]
        ] ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md44", null ],
        [ "Group: math", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md45", [
          [ "math/ fixes", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md46", null ]
        ] ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md47", null ],
        [ "Related notes (outside math/, not changed)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md48", null ],
        [ "Group: lib", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md49", [
          [ "lib group: bug triage", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md50", null ]
        ] ],
        [ "mkh.c renormalization check (requested)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md51", null ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md52", null ],
        [ "Group: maos_skyc", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md53", [
          [ "maos_skyc fixes", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md54", null ]
        ] ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md55", null ],
        [ "Group: cuda", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md56", [
          [ "cuda group: bug review", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md57", null ]
        ] ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md58", null ],
        [ "High-impact items", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md59", null ],
        [ "Group: tools_mex", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md60", [
          [ "tools_mex fixes", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md61", null ]
        ] ],
        [ "Decisions", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md62", null ],
        [ "Build check", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md63", null ]
      ] ],
      [ "Follow-up fixes to the first review", "db/d87/md_docs_2claude__changes__20261120.html#claude_followup", [
        [ "sys/daemonize.c: lock_file_version", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md64", null ],
        [ "tools/scheduler_http.c", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md65", null ],
        [ "Doc block placement", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md66", null ]
      ] ],
      [ "monitor: GTK4 deprecated API replacement", "db/d87/md_docs_2claude__changes__20261120.html#claude_gtk4", [
        [ "Job list (tools/monitor_gui.c)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md67", null ],
        [ "Dialogs (tools/monitor.c)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md68", null ],
        [ "Bugs fixed on the way", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md69", null ]
      ] ],
      [ "Second review of sys, math, lib, maos, skyc, tools, cuda and mex", "db/d87/md_docs_2claude__changes__20261120.html#claude_second_review", [
        [ "Changes that can alter results", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md70", null ],
        [ "sys", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md71", null ],
        [ "math", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md72", null ],
        [ "lib", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md73", null ],
        [ "maos (no result change on valid default configs)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md74", null ],
        [ "skyc and tools", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md75", null ],
        [ "mex", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md76", null ],
        [ "Not changed in this commit", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md77", null ]
      ] ],
      [ "Low-priority lib fixes from the second review", "db/d87/md_docs_2claude__changes__20261120.html#claude_lowprio", null ],
      [ "Open items of the second review", "db/d87/md_docs_2claude__changes__20261120.html#claude_open_items", [
        [ "Changes that can alter results", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md78", null ],
        [ "File format", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md79", null ],
        [ "cuda", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md80", null ],
        [ "maos", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md81", null ],
        [ "skyc", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md82", null ],
        [ "tools", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md83", null ],
        [ "Not changed", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md84", null ]
      ] ],
      [ "Build with a C++ compiler (CC=g++)", "db/d87/md_docs_2claude__changes__20261120.html#claude_cxx", [
        [ "configure.ac", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md85", null ],
        [ "sys/thread.h", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md86", null ],
        [ "sys/common.h: ADDR()", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md87", null ]
      ] ],
      [ "CI: CUDA build in the Linux workflow", "db/d87/md_docs_2claude__changes__20261120.html#claude_ci_cuda", null ],
      [ "Unit tests converted to cmocka", "db/d87/md_docs_2claude__changes__20261120.html#claude_cmocka", [
        [ "Rules for the tests", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md88", null ],
        [ "Infrastructure", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md89", null ],
        [ "Library additions needed by the tests", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md90", null ],
        [ "Programs (test cases)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md91", null ]
      ] ],
      [ "Bugs found by the cmocka tests", "db/d87/md_docs_2claude__changes__20261120.html#claude_testbugs", [
        [ "lib/chol.c", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md92", null ],
        [ "sys/bin.c", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md93", null ],
        [ "Other", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md94", null ],
        [ "Tests", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md95", null ]
      ] ],
      [ "Unit test coverage of sys, math and lib", "db/d87/md_docs_2claude__changes__20261120.html#claude_coverage", [
        [ "New programs (test cases)", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md96", null ],
        [ "Extended programs", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md97", null ],
        [ "Not covered", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md98", null ],
        [ "Library bugs found", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md99", null ]
      ] ],
      [ "Library bugs found by the coverage tests", "db/d87/md_docs_2claude__changes__20261120.html#claude_covbugs", [
        [ "Changes that can alter results", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md100", null ],
        [ "sys", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md101", null ],
        [ "math", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md102", null ],
        [ "lib", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md103", null ],
        [ "Tests", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md104", null ]
      ] ],
      [ "Documentation audit and follow-up fixes", "db/d87/md_docs_2claude__changes__20261120.html#claude_audit", [
        [ "Changes that can alter results", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md105", null ],
        [ "Configuration", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md106", null ],
        [ "configure", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md107", null ],
        [ "Code", "db/d87/md_docs_2claude__changes__20261120.html#autotoc_md108", null ]
      ] ],
      [ "macOS unit test failures", "db/d87/md_docs_2claude__changes__20261120.html#claude_macos", null ],
      [ "UndefinedBehaviorSanitizer", "db/d87/md_docs_2claude__changes__20261120.html#claude_ubsan", null ],
      [ "valgrind", "db/d87/md_docs_2claude__changes__20261120.html#claude_valgrind", null ]
    ] ],
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
"d1/ddc/atm_mk13n50p.html",
"d3/dc7/muv_8h.html#a8a200c66e544ed7d4573c03a5ac55b81",
"d4/d38/parms_8h.html#a1b3586fba1bd11580eac766c9f17c889",
"d4/d38/parms_8h.html#a82aabbaf4f8e34a1e874c07345654fba",
"d4/d38/parms_8h.html#ae6da7f16a611cdd4cfc7d6170af4bbb7",
"d5/ddb/atm_mk25p_g_l_50p_f_a.html",
"d6/d94/maos_2types_8h.html#a63cf2faff093d3cee8e594ad00c56bea",
"d6/d94/maos_2types_8h.html#aecfdde04f1823426a4a6bde295ae6e95",
"d9/df5/skyc_20.html",
"da/dfb/dmath_8h.html#ab48bc9b9c8f923e19dafc70b8fff5de2",
"db/dd5/thread_8h.html#ac484ce8461053d81b280ae84d3d2c110",
"dd/dc1/lmath_8h.html#a57fd9bfee0aa615e6ed2ab9ed0fece94",
"de/d07/cmath_8h.html#ae9103b8a83d73f56aac1fb8650afa6d9",
"de/de7/smath_8h.html#aa1ebd8e1d073034133fd5a178b624698",
"df/def/type_8h.html#a33f7b985bf7933943eed9bb36b4676ed",
"df/def/type_8h.html#ace7f877fd34c3f699d35a6a81c081a26",
"globals_l.html"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';