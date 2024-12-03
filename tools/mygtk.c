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
GtkWidget *tool_button_new(int toggle, const char *name, GdkPixbuf *buf){
	GtkWidget *item=NULL;
	if(name || buf){
		item=toggle?gtk_toggle_button_new():gtk_button_new();
		GtkWidget *image=NULL;
		if(buf){
#if GTK_MAJOR_VERSION<=3
			image=gtk_image_new_from_pixbuf(buf);
#else
			image=gtk_image_new_from_paintable(GDK_PAINTABLE(gdk_texture_new_for_pixbuf(buf)));
#endif
		}else{
			image=gtk_image_new_from_icon_name(name
#if GTK_MAJOR_VERSION<=3
			, GTK_ICON_SIZE_SMALL_TOOLBAR
#endif
			);
		}
#if GTK_MAJOR_VERSION<=3
			gtk_button_set_relief(GTK_BUTTON(item), GTK_RELIEF_NONE);
			gtk_button_set_image(GTK_BUTTON(item), image);
#else
			gtk_button_set_has_frame(GTK_BUTTON(item), FALSE);
			gtk_button_set_child(GTK_BUTTON(item), image);
#endif
	} else{
		item=gtk_vseparator_new();
		gtk_widget_set_opacity(item, 0);
#if GTK_MAJOR_VERSION>=3
		gtk_widget_set_hexpand(item, toggle);
#endif
	}
	return item;
}
//todo: consolidate with drawdaemon_gui.new_tool()
//toolbar is actually an hbox
GtkWidget* new_toolbar_item(GtkWidget *toolbar, GtkWidget *item, int toggle, const char *iconname, 
	GdkPixbuf *iconbuf, const char *cmdname, void(*callback)(GtkWidget*, gpointer data), gpointer data){
	if(!item){
		item=tool_button_new(toggle, iconname, iconbuf);
	}
	if(callback){
		g_signal_connect(item, toggle?"toggled":"clicked", G_CALLBACK(callback), data);
	}
	if(cmdname) gtk_widget_set_tooltip_text(item, cmdname);
	//gtk_widget_set_size_request(item, 16, 16);
#if GTK_MAJOR_VERSION<3
	box_append(GTK_BOX(toolbar), item, FALSE, FALSE, 0);
#else
	if(G_TYPE_CHECK_INSTANCE_TYPE(toolbar, GTK_TYPE_BOX)){
		box_append(GTK_BOX(toolbar), item, FALSE, FALSE, 0);
	}else{
	gtk_header_bar_pack_start(GTK_HEADER_BAR(toolbar), item);
	}
#endif
	return item;
}
