#!/bin/sh
#Test driver used by make check. Runs the automake test-driver, then updates
#the unit test status report (docs/unittest_status.md in the build folder) so that
#it is current even when a test fails.
#Usage: unittest_driver.sh TOP_SRCDIR TOP_BUILDDIR [test-driver options] -- PROGRAM
top_srcdir="$1"
top_builddir="$2"
shift 2
$SHELL "$top_srcdir/m4/test-driver" "$@"
ans=$?
#Serialize report updates when tests run in parallel. mkdir is atomic.
lock=unittest_report.lock
count=0
while ! mkdir "$lock" 2>/dev/null ;do
	count=$((count+1))
	if [ $count -gt 600 ];then #stale lock
		rmdir "$lock" 2>/dev/null
	fi
	sleep 0.1 2>/dev/null || sleep 1
done
$SHELL "$top_srcdir/test/unittest_report.sh" "$top_srcdir" "$top_builddir"
rmdir "$lock"
exit $ans
