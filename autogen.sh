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
if ! which autoreconf; then
	read -p "Required packages are not installed, would like to install them? (Y/n):" ans
	if [ "$ans" = "" -o "$ans" = "Y" -o "$ans" = "y" ]; then
		common="autoconf automake libtool make cmake"
		case "`uname`" in
		Linux)
			if which apt; then
				sudo apt install $common gcc git tar bzip2 libfftw3-dev liblapack3 libcmocka-dev libwebsockets-dev
			elif which dnf; then
				sudo dnf install $common gcc git tar bzip2 fftw-devel lapack libcmocka-devel libwebsockets-devel
			else
				echo "System: Linux. Not sure how to install the packages."
			fi
			;;
		Darwin)
			if which brew; then
				brew install $common
			elif which port; then
				sudo port install $common
			else
				echo "System: macOS. Not sure how to install the packages."
			fi
			;;
		*)
			echo "Unknown system: `uname`"
		esac
	fi
fi
autoreconf -i
