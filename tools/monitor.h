/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_MONITOR_H
#define AOS_MONITOR_H
#include <glib.h>
#include <gtk/gtk.h>
#include <glib/gprintf.h>
#include "../sys/sys.h"
/**
   Information for each process.
*/
typedef struct PROC_T{
    int hid;/*host id. hosts[hid] gives hostname. */
    int pid;
    int done;/*done with it. */
    double frac;
    int oldinfo;
    int iseed_old;
    char *path;
    STATUS_T status;
    GtkWidget *entry_start;
    GtkWidget *entry_pid;
    GtkWidget *entry_path;
    GtkWidget *entry_errlo;
    GtkWidget *entry_errhi;
    GtkWidget *entry_iseed;
    GtkWidget *entry_timing;
    GtkWidget *hbox, *vbox;
    GtkWidget *label;
    GtkWidget *btn;
    GtkWidget *btnimage;
    char time3[80];
    GtkTreeRowReference *row;
    gulong btnhandler;
    struct PROC_T *next;
}PROC_T;

/*void proc_remove(int id,int pid, int flag); */
void refresh(PROC_T *p);
void kill_job(PROC_T *p);
void kill_job_event(GtkWidget *btn, GdkEventButton *event, PROC_T *p);
void notify_user(PROC_T *p);
int scheduler_cmd(int host,int pid, int command);
extern GdkColor blue;
extern GdkColor green;
extern GdkColor red;
extern GdkColor yellow;
extern GdkColor white;
extern GdkColor *bg;
extern GdkColor color_even;
extern GdkColor color_odd;

extern PROC_T **pproc;
extern int *nproc;
extern GtkWidget *notebook;
extern GtkWidget **pages;
extern GdkPixbuf *icon_main;
extern GdkPixbuf *icon_finished;
extern GdkPixbuf *icon_failed;
extern GdkPixbuf *icon_running;
void proc_remove(int id,int pid);
void proc_remove_all(int id);
GtkWidget *new_page(int ihost);
void remove_entry(PROC_T *p);
GtkWidget *monitor_new_entry_progress(void);
GtkWidget *monitor_new_progress(int vertical, int length);
#endif
