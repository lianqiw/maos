which libtoolize>/dev/null && libtoolize --copy --force || glibtoolize --copy --force
aclocal --force -I m4
autoheader --force
automake --add-missing --copy --force-missing
autoconf --force
