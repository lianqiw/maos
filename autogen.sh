#aclocal --force -I m4 &&\
#libtoolize --force --copy&&\
#autoconf --force &&\
#autoheader --force &&\
#automake --add-missing --copy --force-missing
#gtkdocize --copy
#use glibtoolize in mac.

#autoreconf --force --install -I m4
#aclocal --force -I m4
aclocal --force -I m4
autoheader --force
automake --add-missing --copy --force-missing
autoconf --force
