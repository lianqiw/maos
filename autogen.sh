#Autotools includes two Core Packages: GNU Autoconf and GNU Automake
#GNU Autoconf include the following tools
# autoconf:   Creates configure from configure.ac
# autoheader: Creates config.h.in from configure.ac
# autoreconf: Run all tools in the right order
# autoscan:   Scan sources for common portability problems and related macros missing from configure.ac
# autoupdate: Update obsolete macros in configure.ac
# ifnames:    Gather identifiers from all #if/#ifdef/... directives
# autom4te:   The heart of Autoconf. It drives M4 and implements the features used by most of the above tools.

#GNU Automake includes the following tools
# automake:   Create Makefile.in from Makefile.am and configure.ac
# aclocal:    Scan configure.ac for uses of third-party macros and gather definitions in aclocal.m4

#GNU M4 is the real macro processor. 

#Use autoreconf to setup the package initially
autoreconf -i
