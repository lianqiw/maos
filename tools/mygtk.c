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
#include "mygtk.h"

//todo: consolidate with drawdaemon_gui.new_tool()
//toolbar is actually an hbox
void new_toolbar_item(GtkWidget *toolbar, const char *iconname, GdkPixbuf *iconbuf, const char *cmdname, void(*callback)(GtkButton *, gpointer data), int data){
	GtkWidget *item;
	if(cmdname){
#if GTK_MAJOR_VERSION<3
		GtkWidget *image;
		if(iconbuf){
			image=gtk_image_new_from_pixbuf(iconbuf);
		} else{
			image=gtk_image_new_from_icon_name(iconname, GTK_ICON_SIZE_SMALL_TOOLBAR);
		}
		item=gtk_button_new();
		gtk_button_set_image(GTK_BUTTON(item), image);
#else
		item=button_new(iconname);
#if GTK_MAJOR_VERSION<=3
		gtk_button_set_relief(GTK_BUTTON(item), GTK_RELIEF_NONE);
#else
		gtk_button_set_has_frame(GTK_BUTTON(item), FALSE);
#endif
		if(iconbuf){
#if GTK_MAJOR_VERSION<=3
			GtkWidget *image=gtk_image_new_from_pixbuf(iconbuf);
			gtk_button_set_image(GTK_BUTTON(item), image);
#else
			GtkWidget *image=gtk_image_new_from_paintable(GDK_PAINTABLE(gdk_texture_new_for_pixbuf(iconbuf)));
			gtk_button_set_child(GTK_BUTTON(item), image);
#endif
		}
#endif
		if(callback) g_signal_connect(item, "clicked", G_CALLBACK(callback), GINT_TO_POINTER(data));
		if(cmdname) gtk_widget_set_tooltip_text(item, cmdname);
		//gtk_widget_set_size_request(item, 16, 16);
	} else{
		item=gtk_vseparator_new();
	}
	box_append(GTK_BOX(toolbar), item, FALSE, FALSE, 0);
}
