#!/bin/sh
#To assess coverage, do the following
#1. compile maos with --enable-gcov
#2. Run maos
#3. run this script at the root of the building directory
#4. view out/index.html
lcov --capture --directory . --output-file coverage.info || exit
genhtml coverage.info --output-directory out
