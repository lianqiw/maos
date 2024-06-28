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


#ifndef AOS_TOOLS_MYGTK_H
#define AOS_TOOLS_MYGTK_H
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include <cairo-svg.h>
#include <cairo.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <pango/pango.h>
#define GTK_VERSION_AFTER(major,minor) (GTK_MAJOR_VERSION>major || (GTK_MAJOR_VERSION==major && GTK_MINOR_VERSION>=minor))
#define GLIB_VERSION_AFTER(major,minor) (GLIB_MAJOR_VERSION>major || (GLIB_MAJOR_VERSION==major && GLIB_MINOR_VERSION>=minor))

#if GTK_MAJOR_VERSION >2
#define gtk_hbox_new(A,B) gtk_box_new(GTK_ORIENTATION_HORIZONTAL, B)
#define gtk_vbox_new(A,B) gtk_box_new(GTK_ORIENTATION_VERTICAL, B)
#define gtk_vseparator_new() gtk_separator_new(GTK_ORIENTATION_VERTICAL)
#define gtk_hseparator_new() gtk_separator_new(GTK_ORIENTATION_HORIZONTAL)
#endif

#if GTK_VERSION_AFTER(3,2)
#define gtk_hscale_new(A,B...) gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, B)
#define gtk_hscale_new_with_range(A...) gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, A)
#endif
#if GTK_VERSION_AFTER(3,4)
#define grid_attach(A,B,C,D,E,F,G) ({gtk_grid_attach(GTK_GRID(A),B,C,D,E,F); if(G) g_object_set(B, "expand", TRUE, NULL);})
#define gtk_table_new(A,B,C) gtk_grid_new()
#define gtk_table_resize(A...)
#else
#define grid_attach(A,B,C,D,E,F,G) gtk_table_attach(GTK_TABLE(A),B,C,D,E,F, (GtkAttachOptions)G, (GtkAttachOptions)0, 0, 0)
#endif


#if !GLIB_VERSION_AFTER(2,32) //GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
#define g_thread_new(A,B,C) g_thread_create(B, C, 0, NULL);
#endif

#ifndef G_VALUE_INIT
#define G_VALUE_INIT {0,{{0}}}
#endif


#if GTK_MAJOR_VERSION<4 //gtk3 or gtk2
#define box_append(box, child, expand, fill, padding) gtk_box_pack_start(GTK_BOX(box), child, expand, fill, padding)
#define box_prepend(box, child, expand, fill, padding) gtk_box_pack_end(GTK_BOX(box), child, expand, fill, padding)
#define button_new(iconname) gtk_button_new_from_icon_name(iconname, GTK_ICON_SIZE_BUTTON)
#define set_cursor(widget, cursor) gdk_window_set_cursor(gtk_widget_get_window(widget), cursor);
#define toggle_button_set_active(widget, value) gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(widget),value)
#define toggle_button_get_active(widget) gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(widget))
#define scrolled_window_new() gtk_scrolled_window_new(NULL, NULL)
#define window_new() gtk_window_new(GTK_WINDOW_TOPLEVEL)
#define entry_set_text(entry, text) gtk_entry_set_text(GTK_ENTRY(entry), text)
#define cursor_new_from_name(name) gdk_cursor_new_from_name(gdk_display_get_default(), name)
#define check_button_set_active(button, value) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), value)
#define check_button_get_active(button) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))
#if GTK_VERSION_AFTER(3,10)
#define window_destroy(window) gtk_window_close(GTK_WINDOW(window))
#else
#define window_destroy(window) gtk_object_destroy(GTK_OBJECT(window))
#endif
#else //gtk4
#define box_append(box, child, ...) gtk_box_append(GTK_BOX(box), child)
#define box_prepend(box, child, ...) gtk_box_prepend(GTK_BOX(box), child)

#define button_new(iconname)  gtk_button_new_from_icon_name(iconname)
#define set_cursor gtk_widget_set_cursor
#define toggle_button_set_active(widget,value) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget),value)
#define toggle_button_get_active(widget) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget))
#define window_destroy(window) gtk_window_destroy(GTK_WINDOW(window))
#define scrolled_window_new() gtk_scrolled_window_new()
#define window_new() gtk_window_new()
#define entry_set_text(entry, text) gtk_entry_buffer_set_text(gtk_entry_get_buffer(GTK_ENTRY(entry)), text, -1)
#define cursor_new_from_name(name) gdk_cursor_new_from_name(name, NULL)
#define check_button_set_active(button, value) gtk_check_button_set_active(GTK_CHECK_BUTTON(button), value)
#define check_button_get_active(button) gtk_check_button_get_active(GTK_CHECK_BUTTON(button))
#endif

void new_toolbar_item(GtkWidget *toolbar, const char *iconname, GdkPixbuf *iconbuf, const char *cmdname, void(*func)(GtkButton *, gpointer data), int data);
#endif

