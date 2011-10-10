which libtoolize>/dev/null && libtoolize --copy || glibtoolize --copy
aclocal --force -I m4
autoheader --force
automake --add-missing --copy --force-missing
autoconf --force
