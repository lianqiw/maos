which libtoolize>/dev/null && libtoolize --copy --force || glibtoolize --copy --force
aclocal --force -I m4
autoheader --force
automake --add-missing --copy --force-missing
sed 's/-fno-common//g' m4/libtool.m4 > m4/libtool_new.m4
mv m4/libtool_new.m4 m4/libtool.m4
autoconf --force
