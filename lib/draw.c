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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <errno.h>
#include <search.h>
#include "../math/mathdef.h"
#include "draw.h"
#include "cure.h"
int DRAW_ID=0;
int DRAW_DIRECT=0;
int disable_draw=0; /*if 1, draw will be disabled  */
PNEW(lock);
#define MAXDRAW 1024
static int sock_helper=-1;
int listening=0;
int draw_single=0;//1: Only draw active frame. 0: draw all frames.
int draw_changed=0; //1: switched pages
long group=0;
/*If not null, only draw those that match draw_fig and draw_fn*/
/**
   Contains functions for data visualization. 

   2013-02-20 

   There are three ways to enable drawing with drawdaemon.

   1) monitor (remote machine) opens a socket to connect to scheduler to request
   showing progress of a maos run. The scheduler passes the socket to maos. The
   monitor passes the other end of the socket to drawdaemon by fork()+exec(). In
   the end maos has direct connect to drawdaemon using sockets. draw() will
   write to the sockets the information to display.

   2) If DRAW_DIRECT=0, get_drawdaemon() will talk to draw_helper() (a
   fork()'ed process) with pid through a socket pair (AF_UNIX type). The
   draw_helper() will then create another socket_pair, pass one end of it to
   draw(), and the other end to drawdaemon() by fork()+exec().

   3 If DRAW_DIRECT=1, get_drawdaemon() will launch drawdaemon directly.

*/
/* List of list*/
typedef struct list_t{
    char *key;
    struct list_t *next;
    struct list_t *child;//child list
}list_t;
int sock_ndraw=0;//number of displays
int sock_ndraw2=0;//number of valid displays
typedef struct{
    int fd;
    int pause;
    list_t *list;
    char *figfn[2];
}sockinfo_t;
sockinfo_t sock_draws[MAXDRAW];

#define CATCH(A) if(A){\
	warning("stwrite to %d failed with %s\n", sock_draw, strerror(errno)); \
	draw_remove(sock_draw,0);					\
	goto end;							\
    }
#define STWRITESTR(A) CATCH(stwritestr(sock_draw,A))
#define STWRITEINT(A) CATCH(stwriteint(sock_draw,A))
#define STWRITECMDSTR(cmd,str) CATCH(stwriteint(sock_draw,cmd) || stwritestr(sock_draw,str))
#define STWRITE(p,len) CATCH(stwrite(sock_draw,p,len));

/**
   Listen to drawdaemon for update of fig, fn. The hold values are stored in figfn.
*/
static void* listen_drawdaemon(sockinfo_t *sock_data){
    listening=1;
    int sock_draw=sock_data->fd;
    char **figfn=sock_data->figfn;
    //info("draw is listening to drawdaemon at %d\n", sock_draw);
    int cmd;
    while(sock_data->fd!=-1 && !streadint(sock_draw, &cmd)){
	switch(cmd){
	case DRAW_FIGFN:
	    {
		char *fig=0, *fn=0;
		streadstr(sock_draw, &fig);
		streadstr(sock_draw, &fn);
		if(figfn[0] && figfn[1] && (strcmp(figfn[0], fig) || strcmp(figfn[1], fn))){
		    draw_changed=1;
		    //info("draw %d switch to fig=%s, fn=%s\n", sock_draw, fig, fn);
		}
		free(figfn[0]);
		free(figfn[1]);
		figfn[0]=fig;
		figfn[1]=fn;
	    }
	    break;
	case DRAW_PAUSE:
	    sock_data->pause=1;
	    info("draw %d paused\n", sock_draw);
	    break;
	case DRAW_RESUME:
	    sock_data->pause=0;
	    info("draw %d resumed\n", sock_draw);
	    break;
	default:
	    warning("cmd=%d is not understood\n", cmd);
	    break;
	}
    }
    //info("draw stop lisening to drawdaemon at %d\n", sock_draw);
    listening=0;
    return NULL;
}

static int list_search(list_t **head, list_t **node, const char *key, int add){
    list_t *p=0;
    for(p=*head; p; p=p->next){
	if(!strcmp(p->key, key)){
	    break;
	}
    }
    int ans=p?1:0;
    if(add && !p){
	p=mycalloc(1,list_t);
	p->key=strdup(key);
	p->next=*head;
	*head=p;
    }
    if(node) *node=p;
    return ans;
}
static void list_destroy(list_t **head){
    list_t *p;
    for(p=*head; p; p=*head){
	*head=p->next;
	if(p->child){
	    list_destroy(&p->child);
	}
	free(p->key);
	free(p);
    }
}
/**Add fd to list of drawing socks*/
int draw_add(int fd){
    if(fd==-1) return -1;
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	if(sock_draws[ifd].fd<0){//fill a empty slot
	    sock_draws[ifd].fd=fd;
	    sock_ndraw2++;
	    thread_new((thread_fun)listen_drawdaemon, &sock_draws[ifd]);
	    disable_draw=0;
	    return 0;
	}
    }
    if(sock_ndraw<MAXDRAW){
	memset(&sock_draws[sock_ndraw], 0, sizeof(sockinfo_t));
	sock_draws[sock_ndraw].fd=fd;
	thread_new((thread_fun)listen_drawdaemon, &sock_draws[sock_ndraw]);
	sock_ndraw++;
	sock_ndraw2++;
	disable_draw=0;
	return 0;
    }else{
	return -1;
    }
}
static void draw_remove(int fd, int reuse){
    if(fd<0) return;
    if(reuse){
	scheduler_send_socket(fd, DRAW_ID);
    }
    close(fd);
    int found=0;
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	if(sock_draws[ifd].fd==fd){
	    found=1;
	    sock_draws[ifd].fd=-1;
	    list_destroy(&sock_draws[ifd].list);
	    free(sock_draws[ifd].figfn[0]);sock_draws[ifd].figfn[0]=NULL;
	    free(sock_draws[ifd].figfn[1]);sock_draws[ifd].figfn[1]=NULL;
	}
    }
  
    if(found){
	sock_ndraw2--;
    }else{
	dbg("draw_remove: fd=%d is not found\n", fd);
    }
    if(sock_ndraw2<=0){
	disable_draw=1;//do not try to restart drawdaemon
    }
}
static int launch_drawdaemon(){
    int sv2[2];
    /*one end of sv2 will be passed back to draw_helper, the other end of sv2
      will be passed to drawdaemon.*/
    if(!socketpair(AF_UNIX, SOCK_STREAM, 0, sv2)){
	char arg1[20];
	snprintf(arg1, 20, "%d", sv2[1]);
	if(spawn_process("drawdaemon", arg1, NULL)<0){
	    warning("spawn drawdaemon failed\n");
	    close(sv2[0]);
	    sv2[0]=-1;
	}
	close(sv2[1]);
    }else{
	perror("socketpair");
	warning("socket pair failed, cannot launch drawdaemon\n");
	disable_draw=1;
	sv2[0]=-1;
    }
    return sv2[0];
}

/**
   A helper routine that forks in the early stage to launch drawdaemon.  We
   don't use the main routine to launch drawdaemon because it may take very long
   to fork a process that consumes Gig's of memory.
*/
void draw_helper(void){
    if(DRAW_DIRECT){
	info("DRAW_DIRECT=1, skip draw_helper\n");
	return;
    }
    if(sock_helper>-1){
	info("draw_helper is already running.\n");
	return;
    }
    int sv[2];
    if(socketpair(AF_UNIX, SOCK_STREAM, 0, sv)){
	perror("socketpair");
	warning("socketpair failed, disable drawing.\n"); 
	disable_draw=1;
    }else{
	info("draw_helper started\n");
    }
    pid_t pid=fork();
    if(pid<0){
	close(sv[0]); 
	close(sv[1]);
	warning("Fork failed. Return.\n");
    }
    if(pid){/*Parent */
	sock_helper=sv[0];
	close(sv[1]);
    }else{/*Child */
	close(sv[0]);
	int cmd;
	while(!streadint(sv[1], &cmd)){
	    int sd=launch_drawdaemon();
	    if(stwritefd(sv[1], sd)){
		warning("write %d to sv[1] failed\n", sd);
	    }
	}
	close(sv[1]);
	_exit(0);//shall never return.
    }
}
/**
   Open a connection to drawdaemon. sock_draw may be set externally, in which case helper=-1.
*/
static int get_drawdaemon(){
    signal(SIGPIPE, SIG_IGN);
    if(sock_ndraw2){//drawdaemon already connected
	return 0;
    }
    if(disable_draw){
	return -1;
    }
    char *display=getenv("DISPLAY");
    if(display && !strlen(display)){//display is not set, we ask scheduler to open a drawdaemon.
	display=0;
    }
    if(DRAW_ID<=0){
	DRAW_ID=getsid(0);
	if(DRAW_ID<0){
	    DRAW_ID=1;
	}
    }
    int sock=-1;
    //First try reusing existing idle drawdaemon
    if(scheduler_recv_socket(&sock, display?DRAW_ID:0)){
	sock=-1;
    }else{//test whether drawdaemon is still running
	if(stwriteint(sock, DRAW_FINAL)){
	    dbg("received socket=%d is already closed.\n", sock);
	    close(sock);
	    sock=-1;
	}
    }
    if(sock==-1 && display){
	if(DRAW_DIRECT || sock_helper<=-1){//directly fork and launch
	    TIC;tic;
	    sock=launch_drawdaemon();
	    toc("Directly launch drawdaemon");
	}else{//use helper to launch
	    if(stwriteint(sock_helper, DRAW_ID) || streadfd(sock_helper, &sock)){
		sock=-1;
		disable_draw=1;
		close(sock_helper);
		sock_helper=-1;
		warning("Unable to talk to the helper to launch drawdaemon\n");
	    }
	    dbg("launch using sock helper: sock=%d\n", sock);
	}
    }
    if(sock!=-1){
	draw_add(sock);
    }
    if(sock_ndraw2>0){
	return 0;
    }else{
	warning("Failed to open drawdaemon\n");
	return 1;
    }
}

/* 
   Check whether we need to draw current page. True if
   1) not yet exist in the list. it is inserted into the list if add is valid and set as current page.
   2) is current page.
   3) when draw_single is not set.
   4) pause must not be set set.
*/
static int check_figfn(int ifd,  const char *fig, const char *fn, int add){
    if(disable_draw || sock_draws[ifd].fd==-1 || sock_draws[ifd].pause) return 0;
    if(!draw_single) return 1;
    list_t *child=0;
    list_search(&sock_draws[ifd].list, &child, fig, add);
    int found=0;
    if(child){
	found=list_search(&child->child, NULL, fn, add);
    }
    char **figfn=sock_draws[ifd].figfn;
    return (!found || (!mystrcmp(figfn[0], fig) && !mystrcmp(figfn[1], fn)));
}
/**
   Tell drawdaemon that this client will no long use the socket. Send the socket to scheduler for future reuse.
*/
void draw_final(int reuse){
    LOCK(lock);
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	int sock_draw=sock_draws[ifd].fd;
	if(sock_draw==-1) continue;
	STWRITEINT(DRAW_FINAL);
      end:
	draw_remove(sock_draw, reuse);
    }
    UNLOCK(lock);
}

/**
   Check whether what we are drawing is current page of any drawdaemon.
*/
int draw_current(const char *fig, const char *fn){
    if(disable_draw) return 0;
    if(!draw_single) return 1;
    int current=0;
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	int sock_draw=sock_draws[ifd].fd;
	if(sock_draw==-1) continue;
	/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	if(check_figfn(ifd, fig, fn, 0)){
	    current=1;
	}
    }
    return current;
}
/**
   Plot the coordinates ptsx, ptsy using style, and optionally plot ncir circles.
*/
int plot_points(const char *fig,    /**<Category of the figure*/
		long ngroup,        /**<Number of groups to plot*/
		loc_t **loc,        /**<Plot arrays of loc as grid*/
		const dcell *dc,    /**<If loc isempty, use cell to plot curves*/
		const int32_t *style,/**<Style of each point*/
		const real *limit,/**<x min, xmax, ymin and ymax*/
		const char *xylog,  /**<Whether use logscale for x, y*/
		const dmat *cir,    /**<Data for the circles: x, y origin, radius, and color*/
		const char *const*const legend, /**<ngroup number of char**/
		const char *title,  /**<title of the plot*/
		const char *xlabel, /**<x axis label*/
		const char *ylabel, /**<y axis label*/
		const char *format, /**<subcategory of the plot.*/
		...){
    if(disable_draw) return 0;
    format2fn;
    LOCK(lock);
    int ans=0;
    if(!get_drawdaemon()){
	for(int ifd=0; ifd<sock_ndraw; ifd++){
	    /*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	    int sock_draw=sock_draws[ifd].fd;
	    if(!check_figfn(ifd, fig, fn, 1)) continue;
	    STWRITEINT(DRAW_FLOAT); STWRITEINT(sizeof(real));
	    STWRITEINT(DRAW_START);
	    if(loc){/*there are points to plot. */
		for(int ig=0; ig<ngroup; ig++){
		    STWRITEINT(DRAW_POINTS); STWRITEINT(loc[ig]->nloc);
		    STWRITEINT(2);
		    STWRITEINT(1);
		    STWRITE(loc[ig]->locx, sizeof(real)*loc[ig]->nloc);
		    STWRITE(loc[ig]->locy, sizeof(real)*loc[ig]->nloc);
		}
		if(dc){
		    warning("both loc and dc are specified, ignore dc.\n");
		}
	    }else if(dc){
		if(ngroup!=dc->nx*dc->ny){
		    warning("ngroup and dimension of dc mismatch\n");
		    ngroup=dc->nx*dc->ny;
		}
		for(int ig=0; ig<ngroup; ig++){
		    int nx=0, ny=0;
		    real *p=NULL;
		    if(dc->p[ig]){
			nx=dc->p[ig]->nx;
			ny=dc->p[ig]->ny;
			p=dc->p[ig]->p;
		    }
		    STWRITEINT(DRAW_POINTS);
		    STWRITEINT(nx);
		    STWRITEINT(ny);
		    STWRITEINT(0);
		    if(p){
			STWRITE(p, sizeof(real)*dc->p[ig]->nx*dc->p[ig]->ny);
		    }
		}
	    }else{
		error("Invalid Usage\n");
	    }
	    if(style){
		STWRITEINT(DRAW_STYLE);
		STWRITEINT(ngroup);
		STWRITE(style,sizeof(uint32_t)*ngroup);
	    }
	    if(cir){
		if(cir->nx!=4){
		    error("Cir should have 4 rows\n");
		}
		STWRITEINT(DRAW_CIRCLE);
		STWRITEINT(cir->ny);
		STWRITE(cir->p, sizeof(real)*cir->nx*cir->ny);
	    }
	    if(limit){/*xmin,xmax,ymin,ymax */
		STWRITEINT(DRAW_LIMIT);
		STWRITE(limit, sizeof(real)*4);
	    }
	    if(xylog){
		STWRITEINT(DRAW_XYLOG);
		STWRITE(xylog, sizeof(const char)*2);
	    }
	    if(format){
		STWRITEINT(DRAW_NAME);
		STWRITESTR(fn);
	    }
	    if(legend){
		STWRITEINT(DRAW_LEGEND);
		for(int ig=0; ig<ngroup; ig++){
		    STWRITESTR(legend[ig]);
		}
	    }
	    STWRITECMDSTR(DRAW_FIG,fig);
	    STWRITECMDSTR(DRAW_TITLE,title);
	    STWRITECMDSTR(DRAW_XLABEL,xlabel);
	    STWRITECMDSTR(DRAW_YLABEL,ylabel);
	    STWRITEINT(DRAW_END);
	  end:
	    ans=1;
	}
    }
    UNLOCK(lock); 
    return ans;
}

/**
   Draw an image.
*/
//Data structure for imagesc
#define MAXNX 1000 //downsample
typedef float dtype;
typedef struct imagesc_t{
    char *fig; /**<Category of the figure*/
    long nx;   /**<the image is of size nx*ny*/
    long ny;   /**<the image is of size nx*ny*/
    long bsize;/**<bytes of each element*/
    dtype *limit; /**<x min; xmax; ymin and ymax*/
    dtype *zlim;/**< min;max; the data*/
    dtype *p; /**<The image*/
    char *title;  /**<title of the plot*/
    char *xlabel; /**<x axis label*/
    char *ylabel; /**<y axis label*/
    char *fn;
}imagesc_t;
static void* imagesc_do(imagesc_t *data){
    if(disable_draw) return NULL;
    LOCK(lock);
    if(!get_drawdaemon()){
	char *fig=data->fig;
	long nx=data->nx;
	long ny=data->ny;
	dtype *limit=data->limit;
	dtype *zlim=data->zlim;
	dtype *p=data->p;
	char *title=data->title;
	char *xlabel=data->xlabel;
	char *ylabel=data->ylabel;
	char *fn=data->fn;
	for(int ifd=0; ifd<sock_ndraw; ifd++){
	    /*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	    int sock_draw=sock_draws[ifd].fd;
	    if(!check_figfn(ifd, fig, fn, 1)) continue;
	    int32_t header[2];
	    header[0]=nx;
	    header[1]=ny;
	    STWRITEINT(DRAW_FLOAT);
	    STWRITEINT(sizeof(dtype));
	    STWRITEINT(DRAW_START);
	    STWRITEINT(DRAW_DATA);
	    STWRITE(header, sizeof(int32_t)*2);
	    STWRITE(p, sizeof(dtype)*nx*ny);
	    if(zlim){
		STWRITEINT(DRAW_ZLIM);
		STWRITE(zlim,sizeof(dtype)*2);
	    }
	    if(limit){/*xmin,xmax,ymin,ymax */
		STWRITEINT(DRAW_LIMIT);
		STWRITE(limit, sizeof(dtype)*4);
	    }
	    if(fn){
		STWRITECMDSTR(DRAW_NAME,fn);
	    }
	    STWRITECMDSTR(DRAW_FIG,fig);
	    STWRITECMDSTR(DRAW_TITLE,title);
	    STWRITECMDSTR(DRAW_XLABEL,xlabel);
	    STWRITECMDSTR(DRAW_YLABEL,ylabel);
	    STWRITEINT(DRAW_END);
	  end:;
	}
    }
    UNLOCK(lock);
    free(data->fig);
    free(data->limit);
    free(data->zlim);
    free(data->p);
    free(data->title);
    free(data->xlabel);
    free(data->ylabel);
    free(data->fn);
    free(data);
    return NULL;
}
/**
   Draws an image. 

   It uses a separate thread to avoid slowing down the simulation. Skip if socket is busy.
 */
int imagesc(const char *fig, /**<Category of the figure*/
	     long nx,   /**<the image is of size nx*ny*/
	     long ny,   /**<the image is of size nx*ny*/
	     const real *limit, /**<x min, xmax, ymin and ymax*/
	     const real *zlim,/**< min,max, the data*/
	     const real *p, /**<The image*/
	     const char *title,  /**<title of the plot*/
	     const char *xlabel, /**<x axis label*/
	     const char *ylabel, /**<y axis label*/
	     const char *format, /**<subcategory of the plot.*/
	     ...){
    format2fn;
    if(!p || !draw_current(fig, fn)){
	return 0;
    }
    if (draw_single){
	//Skip this drawing if line is busy.
	if((TRYLOCK(lock))){//lock failed
	    return 1;
	}else{
	    UNLOCK(lock);
	}
    }
    //We copy all the data and put the imagesc job into a task
    //The task will free the data after it finishes.
    imagesc_t *data=mycalloc(1, imagesc_t);
#define datastrdup(x) data->x=(x)?strdup(x):0
#define datamemdup(x, size, type)		\
    if(x){					\
	data->x=mymalloc(size, type);		\
	for(long i=0; i<size; i++){		\
	    data->x[i]=(type)x[i];		\
	}					\
    }else{					\
	data->x=0;				\
    }
    datastrdup(fig);
    data->bsize=sizeof(dtype);//always send float.
    datamemdup(limit, 4, dtype);
    datamemdup(zlim, 2, dtype);
    {
	//down sampling and change to float
	int xstep=(nx+MAXNX-1)/MAXNX;
	int ystep=(ny+MAXNX-1)/MAXNX;
	int nx2=(nx)/xstep;
	int ny2=(ny)/ystep;

	data->p=malloc(nx2*ny2*sizeof(dtype));
	for(int iy=0; iy<ny2; iy++){
	    dtype *p2=data->p+iy*nx2;
	    const real *p1=p+iy*ystep*nx;
	    for(int ix=0; ix<nx2; ix++){
		p2[ix]=(dtype)p1[ix*xstep];
	    }
	}
	data->nx=nx2;
	data->ny=ny2;
    }
    datastrdup(title);
    datastrdup(xlabel);
    datastrdup(ylabel);
    data->fn=format?strdup(fn):0;
#undef datastrdup
#undef datamemdup
    thread_new((thread_fun)imagesc_do, data);
    return 1;
}

/**
   Draw the OPD of real and imaginary of complex p defined on nx*ny grid. see imagesc()
*/
int imagesc_cmp_ri(const char *fig, long nx, long ny, const real *limit, const real *zlim,
		    const comp *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;

    real *pr,*pi;
    pr=mymalloc(nx*ny,real);
    pi=mymalloc(nx*ny,real);
    for(int i=0; i<nx*ny; i++){
	pr[i]=creal(p[i]);
	pi[i]=cimag(p[i]);
    }
    imagesc(fig, nx, ny, limit, zlim,pr, title, xlabel, ylabel, "%s real",fn);
    free(pr);
    imagesc(fig, nx, ny, limit, zlim,pi, title, xlabel, ylabel, "%s imag",fn);
    free(pi);
    return 1;
}
/**
   Draw the OPD of abs and phase of complex p defined on nx*ny grid. see imagesc()
*/
int imagesc_cmp_ap(const char *fig, long nx, long ny, const real *limit,const real *zlim,
		    const comp *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    real *pr,*pi;
    int isreal=1;
    pr=mymalloc(nx*ny,real);
    pi=mycalloc(nx*ny,real);
    for(int i=0; i<nx*ny; i++){
	pr[i]=cabs(p[i]);
	if(pr[i]>1.e-10){
	    pi[i]=atan2(cimag(p[i]),creal(p[i]));
	    if(isreal&&fabs(pi[i])>1.e-10) isreal=0;
	}
    }
    imagesc(fig, nx, ny, limit, zlim, pr, title, xlabel, ylabel,
	    "%s abs",fn);
    free(pr);
    if(!isreal){
	imagesc(fig, nx,  ny, limit, zlim, pi,title, xlabel, ylabel,
		"%s phi",fn);
    }
    free(pi);
    return 1;
}
/**
   Draw the OPD of abs of complex p defined on nx*ny grid. see imagesc()
*/
int imagesc_cmp_abs(const char *fig, long nx, long ny, const real *limit,const real *zlim,
		     const comp *p, const char *title, const char *xlabel, const char *ylabel,
		     const char *format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    real *pr;
    pr=mymalloc(nx*ny,real);
    for(int i=0; i<nx*ny; i++){
	pr[i]=cabs(p[i]);
    }
    imagesc(fig, nx, ny, limit, zlim, pr, title, xlabel, ylabel,
	    "%s abs",fn);
    free(pr);
    return 1;
}
#include "../math/mathdef.h"
#include "../math/mathdef.h"
/*
  The following routines applies the imagesc_* functions onto
  dmat,cmat,loc,map,etc, data types. 
   
  The first argument is the type of the plot. Plots with the
  same types are grouped into a single tab.
   
  The last argument is the title of the plot. Plots with the
  same type and title as existing plots will replace the already
  existing plots.
*/
/**
   Mapping the floating point numbers onto screen with scaling similar to matlab
   imagesc.  . see imagesc()
*/

int ddraw(const char *fig, const dmat *A, real *xylim, real *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    return imagesc(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolution value of complex array. see imagesc()
*/
int cdrawabs(const char *fig, const cmat *A, real *xylim, real *zlim,
	      const char *title, const char *xlabel, const char *ylabel,
	      const char *format,...){
  
    format2fn;
    return imagesc_cmp_abs(fig,A->nx,A->ny,xylim,zlim,A->p,title,xlabel,ylabel,"%s abs",fn);
}

/**
   Mapping the real/imaginary part of complex array. see imagesc()
*/
int cdrawri(const char *fig, const cmat *A, real *xylim, real *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...){
    
    format2fn;
    return imagesc_cmp_ri(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolute and phase of complex array. see imagesc()
*/
int cdraw(const char *fig, const cmat *A, real *xylim, real *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    return imagesc_cmp_ap(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   like ddraw, acting on map object. see imagesc()
*/
int drawmap(const char *fig, const map_t *map,  real *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    real limit[4];
    limit[0]=map->ox-map->dx/2;
    limit[1]=map->ox+(map->nx-0.5)*map->dx;
    limit[2]=map->oy-map->dx/2;
    limit[3]=map->oy+(map->ny-0.5)*map->dx;
    imagesc(fig, map->nx, map->ny, limit, zlim, map->p,  title, xlabel, ylabel,"%s",fn);
    return 1;
}
/**
   Plot the loc on the screen. see imagesc()
*/
int drawloc(const char *fig, loc_t *loc, real *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    loc_create_map(loc);
    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    real *opd0=mycalloc(nx*ny,real);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    opd0[ix+iy*nx]=(loc->map->p[(ix+npad)+(iy+npad)*nxm]>0);
	}
    }
    real limit[4];
    limit[0]=loc->map->ox+loc->dx*(npad-1/2);
    limit[1]=limit[0]+loc->dx*nx;
    limit[2]=loc->map->oy+loc->dx*(npad-1/2);
    limit[3]=limit[2]+loc->dx*ny;
    imagesc(fig, nx, ny,limit,zlim,opd0,  title, xlabel, ylabel,"%s",fn);
    free(opd0);
    return 1;
}

/**
   Plot the opd using coordinate loc. see imagesc()
*/
int drawopd(const char *fig, loc_t *loc, const real *opd,  real *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){

    format2fn;
    if(!draw_current(fig, fn)) return 0;
    loc_create_map(loc);
    //This is different from loc_embed. It removes the padding.
    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    dmat *opd0=dnew(nx,ny);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    long ii=loc->map->p[(ix+npad)+(iy+npad)*nxm];
	    if(ii>0){
		opd0->p[ix+iy*nx]=opd[ii-1];
	    }else{
		opd0->p[ix+iy*nx]=NAN;
	    }
	}
    }
    real limit[4];
    limit[0]=loc->map->ox+fabs(loc->dx)*(npad-1/2);
    limit[1]=loc->map->ox+fabs(loc->dx)*(nx+npad-1/2);
    limit[2]=loc->map->oy+fabs(loc->dy)*(npad-1/2);
    limit[3]=loc->map->oy+fabs(loc->dy)*(ny+npad-1/2);
    imagesc(fig,nx,ny, limit,zlim,opd0->p,  title, xlabel, ylabel,"%s",fn);
    dfree(opd0);
    return 1;
}
/**
   Plot gradients using CuReD
*/
int drawgrad(const char *fig, loc_t *saloc, const dmat *grad, real *zlim, int grad2opd,
	      const char *title, const char *xlabel, const char *ylabel,
	      const char* format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    if(grad2opd && grad->nx>8){
	//This is different from loc_embed. It removes the padding.
	dmat *phi=0;
	cure_loc(&phi, grad, saloc);
	real limit[4];
	int npad=saloc->npad;
    	limit[0]=saloc->map->ox+fabs(saloc->dx)*(npad-1/2);
	limit[1]=saloc->map->ox+fabs(saloc->dx)*(phi->nx+npad-1/2);
	limit[2]=saloc->map->oy+fabs(saloc->dy)*(npad-1/2+1);
	limit[3]=saloc->map->oy+fabs(saloc->dy)*(phi->ny+npad-1/2+1);
	//writebin(phi, "phi");
	imagesc(fig,phi->nx,phi->ny, limit, zlim, phi->p,  title, xlabel, ylabel,"%s",fn);
	dfree(phi);

    }else{
	drawopd(fig, saloc, grad->p, zlim, title, xlabel, ylabel, "%s x", fn);
	drawopd(fig, saloc, grad->p+grad->nx/2, zlim, title, xlabel, ylabel, "%s y", fn);
    }
    return 1;
}
/**
   Plot opd*amp with coordinate loc. see imagesc()
*/
int drawopdamp(const char *fig, loc_t *loc, const real *opd, const real *amp, real *zlim,
		const char *title, const char *xlabel, const char *ylabel,
		const char* format,...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    (void)fig;
    loc_create_map(loc);

    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    real ampthres;
    dmaxmin(amp, loc->nloc, &ampthres, 0);
    ampthres*=0.5;
    real *opd0=mycalloc(nx*ny,real);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    long ii=loc->map->p[(ix+npad)+(iy+npad)*nxm]-1;
	    if(ii>-1 && amp[ii]>ampthres){
		opd0[ix+iy*nx]=opd[ii];
	    }else{
		opd0[ix+iy*nx]=NAN;
	    }
	}
    }
    real limit[4];
    limit[0]=loc->map->ox+loc->dx*(npad-1);
    limit[1]=limit[0]+loc->dx*nx;
    limit[2]=loc->map->oy+loc->dx*(npad-1);
    limit[3]=limit[2]+loc->dx*ny;
    imagesc(fig, nx,ny,limit,zlim, opd0,  title, xlabel, ylabel,"%s",fn);
    free(opd0);
    return 1;
}
