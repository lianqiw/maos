/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/

const gchar* all_style=
"progress, trough {\n"
"min-height:4px;"
"min-width: 4px;"
"}\n"
"*{\n"
"padding:1px;" //OK.
"border-radius:1px;" //OK
"border-width:1px;" //OK
"font:12px Sans;"
"}\n"
"notebook{\n"
"background-color:@theme_bg_color;"
"}\n"
"notebook tab{\n" //OK
"border-width: 1px 0px 0px 0px;" //OK
"border-color: #FF0000;"
"border-style: none;"
"border-radius:4px 4px 0px 0px;"  //OK
"padding: 0px 0px 0px 0px;" //OK
"background-color:@selected_bg_color;" //OK
"}\n"
"notebook tab:checked{\n" //OK
"background-color:@theme_bg_color;"
"}\n"
"textview{\n"
"font:14px Sans;"
"border-width:3px;"
"}\n"
".entry.progressbar{\n"
"-GtkEntry-has-frame:1;"
"-GtkEntry-progress-border:0px,0px,0px,0px;"
"-GtkEntry-inner-border:0px,0px,0px,0px;"
"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#0000FF), to (#0000FF));"
"border-width:1px;"
"border-style:none; \n"
"border-color:#000000;"
"}\n"
;
