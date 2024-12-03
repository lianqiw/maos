/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


/*
CSS selector: 
  A B    means B inside A (descendant, includes child, grand-child, etc.).
  A,B    means A or B
  #name: widget name must be prefixed aith a # character
  .class: class names must be prefixed with a period.
In CSS, a margin is the space around an element's border, while padding is the
space between an element's border and the element's content. Put another way,
the margin property controls the space outside an element, and the padding
property controls the space inside an element.

min-height and min-width are affective at controller widget size
//padding is top, right, bottom, left
//border-radius is TL, TR, BR, BL
*/
const gchar* all_style=
//progressbar progress:backdrop is when window is non-active.
"progressbar trough, progress{" //trough is the entire bar. progress is active area(?)
"min-height:15px;"
"min-width: 4px;" 
"border-radius:2px;"
//"background-image: linear-gradient(0deg, yellow, red);"
//"background-color: blue; "
"}"
"*{" //for everything
//"min-height:5px;"
//"min-width: 5px;" //"background-color: rgba(255, 0, 0, 255);"
"margin:0px;" //OK. margin is outside
"padding:1px;" //OK. Padding is blank inside
"border-radius:1px;" //OK
"border-width:1px;" //OK
"font:12px Sans;"
"}"
"entry{"
"min-height:12px;"
"padding: 2px;"
"}"
"checkbutton *{"
"min-height:8px;"
"min-width:8px;"
"}"
"scale *{" //need * to match sub components
"min-height:6px;"
"min-width:6px;"
"}"
"box {"
"padding: 2px;"
"}"
"notebook{\n"
"background-color:@theme_bg_color;"
"}\n"
"notebook tab{\n" //OK
"border-width: 1px;" //OK
//"border-color: #FF0000;"
"border-style: none;"
//"padding: 2px 2px 5px 2px;" 
"background-color:@selected_bg_color;" //OK. selected_bg_color is grayish
"min-height: 8px;"
"min-width:  32px;"
"}\n"
"notebook header.top tab{"
"border-radius: 4px 4px 0px 0px;"  
"}\n"
"notebook header.right tab{"
"border-radius: 0px 4px 4px 0px;" 
"}\n"
"notebook tab:checked, tab:checked *, tab:checked * *{\n" //OK
"background-color:@theme_bg_color;"//theme_bg_color is white
"}\n"
"notebook tab button {"
"padding: 0px;" //OK. Padding is blank inside
"min-height:4px;"
"min-width: 4px;"
"}\n"
"spinbutton, spinbutton *{"
"min-height: 8px;"
"min-width:  8px;"
"}\n"
"textview{\n"
"font:14px Sans;"
"border-width:3px;"
"}\n"
".entry.progressbar{\n"
/*"-GtkEntry-has-frame:1;"
"-GtkEntry-progress-border:0px,0px,0px,0px;"
"-GtkEntry-inner-border:0px,0px,0px,0px;"*/
"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#0000FF), to (#0000FF));"
"border-width:1px;"
"border-style:none;"
"border-color:#000000;"
"}\n"
"headerbar{"
"min-height:8px;"
"border-width: 0px; margin: 0px; padding: 0px;"
"}\n"
"button .toggle {"
"border-width: 0px; margin: 0px; padding: 0px;"
"}\n"
;
