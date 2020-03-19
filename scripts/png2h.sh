#!/bin/sh

for fn in icon*.png;do
    name=${fn/icon-/}
    name=${name/.png/}
    gdk-pixbuf-csource --name=icon_inline_${name} $fn > icon-${name}.h
done
