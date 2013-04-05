/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_TOOLS_MYGTK_H
#define AOS_TOOLS_MYGTK_H
#if GTK_MAJOR_VERSION >2
#define gtk_hbox_new(A,B) gtk_box_new(GTK_ORIENTATION_HORIZONTAL, B)
#define gtk_vbox_new(A,B) gtk_box_new(GTK_ORIENTATION_VERTICAL, B)
#define gtk_vseparator_new() gtk_separator_new(GTK_ORIENTATION_VERTICAL)
#define gtk_hseparator_new() gtk_separator_new(GTK_ORIENTATION_HORIZONTAL)
#endif

#if GTK_MAJOR_VERSION > 2 && GTK_MINOR_VERSION >= 2
#define gtk_hscale_new(A,B...) gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, B)
#define gtk_hscale_new_with_range(A...) gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, A)
#endif
#if GTK_MAJOR_VERSION >=3 && GTK_MINOR_VERSION >= 4
#define grid_attach(A,B,C,D,E,F,G) ({gtk_grid_attach(GTK_GRID(A),B,C,D,E,F); if(G) g_object_set(B, "expand", TRUE, NULL);})
#define gtk_table_new(A,B,C) gtk_grid_new()
#define gtk_table_resize(A...)
#else
#define grid_attach(A,B,C,D,E,F,G) gtk_table_attach(GTK_TABLE(A),B,C,D,E,F, (GtkAttachOptions)G, (GtkAttachOptions)0, 0, 0)
#endif


#if GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
#define g_thread_new(A,B,C) g_thread_create(B, C, 0, NULL);
#endif

#ifndef G_VALUE_INIT
#define G_VALUE_INIT {0,{{0}}}
#endif

#endif

