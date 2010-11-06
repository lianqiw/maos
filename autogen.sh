#aclocal --force -I m4 &&\
#libtoolize --force --copy&&\
#autoconf --force &&\
#autoheader --force &&\
#automake --add-missing --copy --force-missing
#gtkdocize --copy
#use glibtoolize in mac.
autoreconf --force --install -I m4
