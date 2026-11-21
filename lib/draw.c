/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <errno.h>
#include <search.h>

#include "draw.h"
#include "cure.h"
PNEW(lock);//protect the list and TCP socket from race condition
#define MAXDRAW 1024


//Every command is prefixed with DRAW_ENTRY with payload length (not including the command). 
//This helps client to handle unknown command as well as websocket proxy to forward the data.
#define WRITECMD(sock, cmd, npayload) (stwriteint(sock, DRAW_ENTRY)||stwriteint(sock, npayload)||stwriteint(sock, cmd)) 
#define WRITECMDSTR(sock, cmd, str) (WRITECMD(sock, cmd, sizeof(int)+strlen(str)+1)||stwritestr(sock, str))
#define WRITECMDARR(sock, cmd, p, len) (WRITECMD(sock, cmd, len)||stwrite(sock, p, len))
#define WRITECMDINT(sock, cmd, val) (WRITECMD(sock, cmd, sizeof(int))||stwriteint(sock, val))
/**Client identification, usually PID. Same draw_id reuses drawdaemon.*/
int draw_id=DRAW_ID_MAOS;
/**Directly launch drawdaemon without forking a draw_helper*/
int draw_direct=0;
/**if set, draw will be disabled*/
int draw_disabled=0;
/**if ==1, only draw active frame and skip if line is busy. if ==-1, disable draw_single altogether. otherwise draw all frames.*/
int draw_single=0;
/**socket of draw_helper*/
static int sock_helper=-1;

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

   2) If draw_direct=0, get_drawdaemon() will talk to draw_helper() (a
   fork()'ed process) with pid through a socket pair (AF_UNIX type). The
   draw_helper() will then create another socketpair, pass one end of it to
   draw(), and the other end to drawdaemon() by fork()+exec().

   3 If draw_direct=1, get_drawdaemon() will launch drawdaemon directly.

*/

/**
   List of list
*/
typedef struct list_t{
	char* key;
	struct list_t* next;
	struct list_t* child;//child list
}list_t;

typedef struct sockinfo_t{
	int fd;//socket;
	int pause;//do not draw anything
	int draw_single;//only draw active figure
	list_t* list; //list of existing figures
	char* figfn[2];//active figure group and name
    pthread_t thread;//thread that is listening
	struct sockinfo_t* next;//next dradaemon
}sockinfo_t;
sockinfo_t *sock_draws=NULL;
int use_udp=0;
/**
   Write an int to a buffer stream.

   @param[in] fbuf Buffer stream.
   @param[in] A Value to write.
   @return 0 on success, -1 on failure.
*/
static inline int fwriteint(FILE* fbuf, int A){
	if(fwrite(&A, sizeof(int), 1, fbuf)<1) return -1;
	return 0;
}
/**
   write str length and str itself (with tailing NULL)

   @param[in] fbuf Buffer stream.
   @param[in] str String to write. NULL is written as an empty string.
   @return 0 on success, -1 on failure.
*/
static inline int fwritestr(FILE* fbuf, const char* str){
	if(!str) str="";
	uint32_t nlen=(uint32_t)(strlen(str)+1);
	if(fwriteint(fbuf, nlen) || fwrite(str, 1, nlen, fbuf)<nlen) return -1;
	return 0;
}

#define CATCH(A) if(A) {ans=1; warning("fwrite failed\n"); goto end2;}
#define FWRITEARR(p,len) CATCH(fwrite(p,1,len,fbuf)<len)
#define FWRITESTR(str) CATCH(fwritestr(fbuf, str)) //will write 1 byte if str is NULL
#define FWRITEINT(A) CATCH(fwriteint(fbuf,A))
#define FWRITECMDSTR(cmd,str) if(str){FWRITEINT(cmd); FWRITESTR(str);}
#define FWRITECMDARR(cmd,p,len) {FWRITEINT(cmd); if(len>0) FWRITEARR(p,len);}
#define FWRITECMDINT(cmd,val) {FWRITEINT(cmd); FWRITEINT(val);}
//Plot an empty page
#define BUF_POSTPROC\
	fclose(fbuf);\
	if(ans){\
		warning("write failed:%d\n", ans);\
		free_default(buf); bufsize=0; buf=0;\
	} else{\
		int* bufp=(int*)(buf);\
		/*bufp[0] is DRAW_ENTRY*/\
		bufp[1]=(int)bufsize-3*sizeof(int);\
	}

int plot_empty(sockinfo_t *ps, const char *fig, const char *fn);
/**
   Listen to drawdaemon for update of fig, fn. The hold values are stored in figfn.
   Also handles pause, resume, and draw_single toggling. Runs in a separate thread.

   @param[in,out] psock_data Drawdaemon connection; figfn, pause, and draw_single are updated.
   @return NULL.
*/
static void* listen_drawdaemon(void *psock_data){
	sockinfo_t *sock_data=(sockinfo_t *)psock_data;
	int fd=sock_data->fd;
	char** figfn=sock_data->figfn;
	//info("draw is listening to drawdaemon at %d\n", sock_draw);

	int cmd;
	int nlen=0;

	if(!socket_recv_timeout(fd, 0)){
        //make sure it blocks when no data is readable
        dbg("started listening at %d\n", fd);
    }else{
        dbg("set timeout=0 failed for %d\n", fd);
    }
	while(sock_data->fd!=-1){
		int ans=streadint(fd, &cmd);
		if(ans<0){
			warning("streadint %d failed with %d. Stop listening.\n", fd, ans);
			break;
		}
		if(cmd==DRAW_ENTRY){//every message in new format start with DRAW_ENTRY.
			streadint(fd, &nlen);
			streadint(fd, &cmd);
		}
		switch(cmd){
		case DRAW_FIGFN:
		{
			char* fig=0, * fn=0;
			streadstr(fd, &fig);
			streadstr(fd, &fn);
			if(fig&&fn){
				if(figfn[0]&&figfn[1]&&(strcmp(figfn[0], fig)||strcmp(figfn[1], fn))){
					if(sock_data->draw_single==1){
						dbg2("draw %d switch to fig=%s, fn=%s\n", fd, fig, fn);
					}
				}
				free(figfn[0]);
				free(figfn[1]);
				figfn[0]=fig;
				figfn[1]=fn;
			}
		}
			break;
		case DRAW_PAUSE:
			sock_data->pause=1;
			dbg2("draw %d paused\n", fd);
			break;
		case DRAW_RESUME:
			sock_data->pause=0;
			dbg2("draw %d resumed\n", fd);
			break;
		case DRAW_SINGLE:
			sock_data->draw_single=!sock_data->draw_single;
			info("draw %d draw_single change to %d\n", fd, sock_data->draw_single);
			break;
		default:
			if(nlen){
				void* p=malloc(nlen);
				stread(fd, p, nlen);
				free(p);
			}
			dbg("Unknown cmd: %d with size %d from socket %d. Close connection.\n", cmd, nlen, fd);
			//do not close(fd) here: sockinfo_close() closes p->fd later and the number may have been reused.
			shutdown(fd, SHUT_RDWR);
			break;
		}
		cmd=0;
	}
	//dbg_time("Stopped listening to drawdaemon at %d, errno=%d, %s\n", fd, errno, strerror(errno));
	return NULL;
}

/**
   Search for key in a linked list, optionally adding it to the head if not found.

   @param[in,out] head Head of the list. Updated if a node is added.
   @param[out] node If not NULL, set to the found or added node (NULL if not found and not added).
   @param[in] key Key to search for.
   @param[in] add Whether to add the key if not found.
   @return 1 if key is found, 0 otherwise (even if it is added).
*/
static int list_search(list_t** head, list_t** node, const char* key, int add){
	list_t* p=0;
    int ans=0;
	for(p=*head; p; p=p->next){
		if(!strcmp(p->key, key)){
            ans=1;
			break;
		}
	}
	if(add&&!ans){
		p=mycalloc(1, list_t);
		p->key=strdup(key);
		p->next=*head;
		*head=p;
	}
	if(node) *node=p;
	return ans;
}
/**
   Free a linked list and its children recursively.

   @param[in,out] head Head of the list. Set to NULL.
*/
static void list_destroy(list_t** head){
	list_t* p;
	for(p=*head; p; p=*head){
		*head=p->next;
		if(p->child){
			list_destroy(&p->child);
		}
		free(p->key);
		free(p);
	}
}
/**
   Add fd to list of drawing socks. lock is already locked. Sends client
   information to drawdaemon and starts a thread listening to it.

   @param[in] sock Socket connected to drawdaemon. Closed if writing to it fails.
   @return 0 on success, -1 on failure.
*/
int draw_add(int sock){
	//Check that the drawdaemon is live.
	if(sock==-1) return -1;
	draw_disabled=0;
	char *buf=0;
	size_t bufsize=0;
	FILE *fbuf=open_memstream(&buf,&bufsize);
	int ans=0;
	FWRITECMDINT(DRAW_ENTRY, 0);
	FWRITECMDINT(DRAW_PID, getpid());
	if(EXENAME){
		FWRITECMDSTR(DRAW_EXENAME, EXENAME);
	}
	if((DIRSTART||DIROUT)){
		FWRITECMDSTR(DRAW_PATH, DIROUT?DIROUT:DIRSTART);
	}
end2:
	BUF_POSTPROC;
	if(bufsize && !ans && stwrite(sock,buf,bufsize)){
		info("write to %d failed: %s\n",sock,strerror(errno));
		close(sock);
		ans=-1;
	}
	free_default(buf);
	if(ans==-1) return -1;
	sockinfo_t *p=NULL;
	for(p=sock_draws; p; p=p->next){
		if(p->fd==-1){
			break;
		}
	}
	if(!p){
		p=mycalloc(1, sockinfo_t);
		p->next=sock_draws;
		sock_draws=p;
	}
	p->fd=sock;
	if(draw_single!=-1){
		p->draw_single=1;
	}
	pthread_create(&p->thread, NULL, listen_drawdaemon, p);
	return 0;
}
/**
   Close a drawdaemon connection and stop the listening thread. The struct itself is kept.

   @param[in,out] p Drawdaemon connection. fd is set to -1 and lists are freed.
   @param[in] reuse Whether to send the socket to the scheduler for reuse.
*/
static void sockinfo_close(sockinfo_t *p, int reuse){
	if(reuse && p->fd!=-1){
		if(!WRITECMD(p->fd, DRAW_FINAL, 0)){
			dbg("send %d back to scheduler for reuse\n", p->fd);
			scheduler_socket(1, &p->fd, draw_id);
		}
	}
	if(p->thread){
		void *ans;
		if(pthread_cancel(p->thread) || pthread_join(p->thread, &ans)){
			dbg("Unable to cancel or join thread\n");
		}
		p->thread=0;
	}
	if(p->fd!=-1){
		close(p->fd);
		p->fd=-1;
	}
	list_destroy(&p->list);
	FREE(p->figfn[0]);
	FREE(p->figfn[1]);
}

/**
   Launch drawdaemon with one end of a socket pair. Sets draw_disabled if the socket pair cannot be created.

   @return The other end of the socket pair, or -1 on failure.
*/
static int launch_drawdaemon(){
	int sv2[2];//socket pair for communication.
	/*one end of sv2 will be passed back to call, the other end of sv2 will be passed to drawdaemon.*/
	if(!socketpair(AF_UNIX, SOCK_STREAM, 0, sv2)){
		char arg1[20];
		snprintf(arg1, 20, "%d", sv2[1]);
		const char *args[3]={"drawdaemon", arg1, NULL};
		
		if(spawn_process("drawdaemon", args, NULL)<0){
			warning("spawn drawdaemon failed\n");
			close(sv2[0]);
			sv2[0]=-1;
		}
		close(sv2[1]);
	} else{
		perror("socketpair");
		warning("socket pair failed, cannot launch drawdaemon, set draw_disabled=1\n");
		draw_disabled=1;
		sv2[0]=-1;
	}
	return sv2[0];
}

/**
   A helper routine can be called in the early stage to launch drawdaemon. Deprecated.
*/
void draw_helper(void){
	if(draw_direct){
		info("draw_direct=1, skip draw_helper\n");
		return;
	}
	if(sock_helper>-1){
		info("draw_helper is already running.\n");
		return;
	}
	int sv[2];
	if(socketpair(AF_UNIX, SOCK_STREAM, 0, sv)){
		perror("socketpair");
		warning("socketpair failed, disable draw_helper.\n");
		draw_direct=1;
		return;
	} else{
		info("draw_helper started\n");
	}
	pid_t pid=fork();
	if(pid<0){
		close(sv[0]);
		close(sv[1]);
		warning("Fork failed. Return.\n");
		return;
	}
	if(pid){//Parent 
		sock_helper=sv[0];
		close(sv[1]);
	} else{//Child 
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
   Tries to reuse an idle drawdaemon via the scheduler, then to launch one. Sets draw_disabled on failure.

   @return 0 if a drawdaemon is connected, nonzero otherwise.
*/
static int get_drawdaemon(){
	signal(SIGPIPE, SIG_IGN);
	if(sock_draws){//drawdaemon already connected
		return 0;
	}
	if(draw_disabled){
		return -1;
	}
	int DRAW_NOREUSE=0;
	READ_ENV_INT(DRAW_NOREUSE, 0, 1);//if ==1, do not reuse previous drawdaemon
	if(DRAW_NOREUSE){
		draw_id=getppid()%10000;
	}
#if __APPLE__
	const char* display=":0";//always available
#else
	const char* display=getenv("DISPLAY"); 
	if(display&&!strlen(display)) display=NULL;
	if(!display){
		display=getenv("WAYLAND_DISPLAY"); 
		if(display&&!strlen(display)) display=NULL;
	}
#endif
	LOCK(lock);
	int sock=-1;
	//First try reusing existing idle drawdaemon with the same id.
	while(sock==-1 && !scheduler_socket(-1, &sock, draw_id)){
		//test whether received drawdaemon is still running
		if(WRITECMD(sock, DRAW_INIT, 0)){
			dbg("received socket=%d is already closed.\n", sock);
			close(sock);
			sock=-1;
		}else{
			dbg("received socket=%d is ok.\n", sock);
		}
	}
	if(sock==-1	&& !scheduler_socket(0, &sock, draw_id)){
		dbg("launch using scheduler: sock=%d\n", sock);
	}
	if(sock==-1 && display){//try launch locally
		if(draw_direct||sock_helper<=-1){//directly fork and launch
			TIC;tic;
			sock=launch_drawdaemon();
			toc2("Directly launch drawdaemon");
		} else{//use helper to launch
			if(stwriteint(sock_helper, draw_id)||streadfd(sock_helper, &sock)){
				sock=-1;
				draw_disabled=1;
				close(sock_helper);
				sock_helper=-1;
				warning("Unable to talk to the helper to launch drawdaemon. set draw_disabled=1.\n");
			}else{
				dbg("launch using sock helper: sock=%d\n", sock);
			}
		}
	}
	//prevent hang. too small timeout will prevent large data from passing through.
	if(sock==-1 || draw_add(sock) || socket_send_timeout(sock, 60)){
		draw_disabled=1;
		sock=-1;
		warning("Unable to open drawdaemon. Disable drawing. set draw_disabled=1\n");
	}
	UNLOCK(lock);
	return sock==-1;
}
/**
   Tell drawdaemon that this client will no long use the socket. Send the socket to scheduler for future reuse.

   @param[in] reuse Whether to send the sockets to the scheduler for reuse.
*/
void draw_final(int reuse){
	//called from other threads, need to lock
	//may be failed if it is invoked by signal handler while draw is in progress
	if(!sock_draws) return;
	if(!TRYLOCK(lock)){
		for(sockinfo_t **curr=&sock_draws; *curr;){
			sockinfo_t *p=*curr;
			*curr=p->next;
			sockinfo_close(p, reuse);
			FREE(p);
		}
		if(!sock_draws){
			dbg("draw_final. set draw_disabled=1\n");
			draw_disabled=1;//do not try to restart drawdaemon
		}
		UNLOCK(lock);
	}
}

/**
   Check whether we need to draw current page. True if
   1) is current page.
   2) when draw_single is not set.
   3) pause must not be set set.
   If the page does not yet exist in the list, it is inserted and an empty page is sent instead.

   @param[in,out] ps Drawdaemon connection. Its list of pages is updated.
   @param[in] fig Category of the figure.
   @param[in] fn Name of the figure. NULL to check fig only.
   @return 0: do not draw. 1: always draw. 2: draw if lock success.
*/
static int check_figfn(sockinfo_t *ps, const char* fig, const char* fn){
	if(draw_disabled||ps->fd==-1||ps->pause) return 0;
	if(!(draw_single==1 && ps->draw_single==1)) return 1;
	LOCK(lock);
	list_t* child=0;
	list_search(&ps->list, &child, fig, 1);
	int found=0;
	if(child){
		if(fn){
			found=list_search(&child->child, NULL, fn, 1);
		} else{//don't check fn if not set
			found=1;
		}
	}
	int ans=0;//default is false
	if(!found){//page is not found
		//ans=3; //draw without test lock
		plot_empty(ps, fig, fn);
	} else{//page is found
		char** figfn=ps->figfn;
		if(!mystrcmp(figfn[0], fig)){
			if(!fn){
				ans=1;//draw without test lock
			} else if(!mystrcmp(figfn[1], fn)){
				ans=2;//draw if lock success
			}
		}
	}
	UNLOCK(lock);
	return ans;
}

/**
   Check whether what we are drawing is current page of any drawdaemon.

   @param[in] fig Category of the figure.
   @param[in] fn Name of the figure. NULL to check fig only.
   @return 0: do not draw. 1: draw. 2: current page and draw is not busy.
*/
int draw_current(const char* fig, const char* fn){
    int current=0;
	if(draw_disabled){
        current=0;
    }else if(draw_single!=1){
        current=1;
    }else{
		for(sockinfo_t *ps=sock_draws; ps; ps=ps->next){
			if(ps->fd==-1) continue;
            if((current=check_figfn(ps, fig, fn))){
                break;
            }
        }
        if(current==2){//check whether draw is busy
			static int nskip=0;
			static int nplot=0;
            if(TRYLOCK(lock)!=0){//trylock failed
                current=0;
				nskip++;
                if(nplot) dbg2("Skip after plot %d times\n", nplot);
				nplot=0;
            } else{//trylock success
				nplot++;
				if(nskip) dbg2("Resume after skip %d times\n", nskip);
				nskip=0;
                UNLOCK(lock);
            }
        }
    }
	return current;
}
/**
   Same as draw_current() with fn specified by printf-style format.

   @param[in] fig Category of the figure.
   @param[in] format printf style format of the figure name, followed by its arguments.
   @return See draw_current().
*/
int draw_current_format(const char *fig, const char *format,...){
	format2fn;
	return draw_current(fig, fn);
}
/**
   Send an empty page with category fig and name fn to the drawdaemon
   connected via ps so that the page appears in the GUI.

   @param[in,out] ps Drawdaemon connection. Closed on write failure.
   @param[in] fig Category of the figure.
   @param[in] fn Name of the figure.
   @return 0 on success, -1 on write failure (the socket is closed), 1 if serialization fails.
*/
int plot_empty(sockinfo_t *ps, const char *fig, const char *fn){
	if(ps->fd==-1) return -1;
	char *buf=0;
	int ans=0;
	size_t bufsize=0;
	FILE *fbuf=open_memstream(&buf,&bufsize);
	FWRITECMDINT(DRAW_ENTRY, 0);
	FWRITEINT(DRAW_START);
	FWRITECMDINT(DRAW_FLOAT,sizeof(real));
	FWRITECMDSTR(DRAW_FIG,fig);
	FWRITECMDSTR(DRAW_NAME,fn);
	FWRITEINT(DRAW_END);
end2:
	BUF_POSTPROC;
	//do not use send_buf here 
	if(bufsize && !ans && stwrite(ps->fd,buf,bufsize)){
		info("write to %d failed: %s\n",ps->fd,strerror(errno));
		ans=-1;
		sockinfo_close(ps, 0); 
	}
	free_default(buf);
	return ans;
}
/**
   Send the serialized plot data in buf (of bufsize bytes) to all connected
   drawdaemons that need page fig/fn (always send if always is set).

   @param[in] fig Category of the figure.
   @param[in] fn Name of the figure.
   @param[in] buf Serialized plot data.
   @param[in] bufsize Size of buf in bytes.
   @param[in] always Whether to send regardless of whether the page is current.
   @return 0 on success, -1 if writing to any drawdaemon failed.
*/
int send_buf(const char *fig, const char *fn, char *buf, size_t bufsize, int always){
	int ans=0;
	for(sockinfo_t *ps=sock_draws; ps; ps=ps->next){
		/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
		int fd=ps->fd;
		if(fd==-1) continue;
		int needed_i=always?1:check_figfn(ps, fig, fn);
		if(needed_i==2){
			if(TRYLOCK(lock)){//line busy
				needed_i=0;
			}
		}else if(needed_i){
			LOCK(lock);
		}
		if(needed_i){
			block_signal(1);//avoid kill/term signal while holding mutex
			if(use_udp){
				error("To be implemented\n");//use sendmmsg with GSO is fastest.
			} else{
				if(stwrite(fd, buf, bufsize)){
					info("write to %d failed: %s\n", fd, strerror(errno));
					ans=-1;
					sockinfo_close(ps, 0);//keep struct to avoice race condition
				}
			}
			block_signal(0);
			UNLOCK(lock);
		}
	}
	return ans;
}
typedef float dtype;
/**
   Plot an image, points, or curves as specified in opts, and optionally circles.
   Connects to drawdaemon if not yet.

   @param[in] fig Category of the figure
   @param[in] opts Additional options: image, points (loc), curves (dc or dd), circles, limits, style, legend, etc.
   @param[in] title title of the plot
   @param[in] xlabel x axis label
   @param[in] ylabel y axis label
   @param[in] format printf style format of the subcategory (name) of the plot, followed by its arguments.
   @return 0 if succeeded or nothing to draw, -1 if sending failed, 1 if serialization failed.
*/
int draw(const char* fig, const plot_opts opts, const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	if(draw_disabled) return 0;
	format2fn;

	int ans=0;
	if(!get_drawdaemon()){
		int needed=opts.always?1:0;
		for(sockinfo_t *ps=sock_draws; needed!=1 && ps; ps=ps->next){
			/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
			int needed_i=check_figfn(ps, fig, fn);
			if(needed_i){
				//dbg("%s %s: %d\n", fig, fn, needed_i);
				if(!needed||needed==3){
					needed=needed_i;
				}
			}
		}
		char* buf=0;
		size_t bufsize=0;
		if(needed){
			//When using UDP, we need to serialize the data instead of writing like a FILE
			//Use open_memstream to serialize and f_openmem for de-serialize
			//To be able to handle the data in forward/backward compatible way,
			//each segment is composed of 1) the length of bytes (int), 2) command name (int), 3) payload
			FILE* fbuf=open_memstream(&buf, &bufsize);
			FWRITECMDINT(DRAW_ENTRY, 0);
			FWRITEINT(DRAW_START);
			FWRITECMDSTR(DRAW_FIG, fig);
			FWRITECMDSTR(DRAW_NAME, fn);
			FWRITECMDINT(DRAW_FLOAT, sizeof(real));
			if(opts.image){
				if(sizeof(real)!=sizeof(dtype)) FWRITECMDINT(DRAW_FLOAT, sizeof(dtype));
				//FWRITEARR(opts.image->p, nlen2);
#define MAXNX 1024 //downsample if bigger.
				const dmat *p=opts.image;
				long xstep=(NX(p)+MAXNX-1)/MAXNX;
				long ystep=(NY(p)+MAXNX-1)/MAXNX;
				int nx2=(int)(NX(p)/xstep);
				int ny2=(int)(NY(p)/ystep);
				dtype *tmp=mymalloc(nx2*ny2, dtype);
				for(int iy=0; iy<ny2; iy++){
					dtype *pout=tmp+iy*nx2;
					if((opts.image->id&M_REAL)==M_REAL){
						const dmat *pd=dmat_cast(opts.image);
						const real *pin=P(pd)+iy*ystep*NX(pd);
						for(int ix=0; ix<nx2; ix++){
							pout[ix]=(dtype)pin[ix*xstep];
						}
					}else if((opts.image->id&M_FLT)==M_FLT){
						const smat *pd=smat_cast(opts.image);
						const float *pin=P(pd)+iy*ystep*NX(pd);
						for(int ix=0; ix<nx2; ix++){
							pout[ix]=(dtype)pin[ix*xstep];
						}
					}else if((opts.image->id&M_COMP)==M_COMP){//convert complex numbers to real.
#define COMP_TO_REAL(ctype)\
switch(ctype){\
	case 0:\
		for(int ix=0; ix<nx2; ix++){\
			pout[ix]=(dtype)cabs(pin[ix*xstep]);\
		}\
		break;\
	case 1:\
		for(int ix=0; ix<nx2; ix++){\
			pout[ix]=(dtype)atan2(cimag(pin[ix*xstep]), creal(pin[ix*xstep]));\
		}\
		break;\
	case 2:\
		for(int ix=0; ix<nx2; ix++){\
			pout[ix]=(dtype)creal(pin[ix*xstep]);\
		}\
		break;\
	case 3:\
		for(int ix=0; ix<nx2; ix++){\
			pout[ix]=(dtype)cimag(pin[ix*xstep]);\
		}\
		break;\
	}	
						const cmat *pc=cmat_cast(opts.image);
						const comp *pin=P(pc)+iy*ystep*NX(pc);
						COMP_TO_REAL(opts.ctype);
					}else if((opts.image->id&M_ZMP)==M_ZMP){
						const zmat *pc=zmat_cast(opts.image);
						const fcomplex *pin=P(pc)+iy*ystep*NX(pc);
						COMP_TO_REAL(opts.ctype);
					}
				}
				int32_t header[2];
				header[0]=nx2;
				header[1]=ny2;
				size_t nlen1=sizeof(header);
				size_t nlen2=sizeof(dtype)*nx2*ny2;
				FWRITEINT(DRAW_DATA);
				FWRITEARR(header, nlen1);
				FWRITEARR(tmp, nlen2);
				free(tmp);
				if(sizeof(real)!=sizeof(dtype)) FWRITECMDINT(DRAW_FLOAT, sizeof(real));
				if(opts.zlim[0]!=opts.zlim[1]){
					FWRITECMDARR(DRAW_ZLIM, opts.zlim, sizeof(real)*2);
				}
				if(opts.zlog){
					FWRITECMDINT(DRAW_ZLOG, opts.zlog);
				}
			}else{//points or lines
				const int ngroup=(opts.dc && (opts.ngroup>PN(opts.dc)||opts.ngroup==0))?(int)PN(opts.dc):(opts.ngroup?opts.ngroup:1);
				if(opts.loc){/*there are points to plot. */
					for(int ig=0; ig<ngroup; ig++){
						long nlen=opts.loc[ig]->nloc;
						if(opts.maxlen && opts.maxlen<nlen) nlen=opts.maxlen;
						FWRITEINT(DRAW_POINTS);
						FWRITEINT((int)nlen);
						FWRITEINT(2);
						FWRITEINT(1);
						FWRITEARR(opts.loc[ig]->locx, sizeof(real)*nlen);
						FWRITEARR(opts.loc[ig]->locy, sizeof(real)*nlen);
					}
					if(opts.dc){
						warning("both loc and dc are specified, ignore dc.\n");
					}
				} else if(opts.dc || opts.dd){
					for(int ig=0; ig<ngroup; ig++){
						dmat* p=opts.dc?P(opts.dc, ig):opts.dd[ig];
						long nx=NX(p);
						long ny=NY(p);
						if(nx==1 && ny>1){
							nx=ny;
							ny=1;
						}
						if(opts.maxlen&&opts.maxlen<nx) nx=opts.maxlen;
						FWRITEINT(DRAW_POINTS);
						FWRITEINT((int)nx);//number of points
						FWRITEINT((int)ny);//number of numbers per point. 1 or 2.
						FWRITEINT(0);//square plot or not
						if(nx && ny){
							FWRITEARR(P(p), nx*ny*sizeof(real));
						}
					}
				} else {
					warning("Empty plot.\n");
				}
				if(opts.dx){
					FWRITEINT(DRAW_DX);
					FWRITEINT(1);
					FWRITEARR(&opts.dx, sizeof(real));
				}else if(opts.dxs){
					FWRITEINT(DRAW_DX);
					FWRITEINT((int)PN(opts.dxs));
					FWRITEARR(P(opts.dxs), sizeof(real)*PN(opts.dxs))
				}
				if(opts.xylog){
					FWRITECMDARR(DRAW_XYLOG, opts.xylog, sizeof(char)*2);
				}
				if(opts.style){
					FWRITEINT(DRAW_STYLE);
					FWRITEINT(ngroup);
					FWRITEARR(opts.style, sizeof(uint32_t)*ngroup);
				}
				if(opts.legend){
					FWRITEINT(DRAW_LEGEND);
					for(int ig=0; ig<ngroup; ig++){
						FWRITESTR(opts.legend[ig]);
					}
				}
			}
			//common
			if(opts.cir){
				if(NX(opts.cir)!=4){
					error("Cir should have 4 rows\n");
				}
				FWRITEINT(DRAW_CIRCLE);
				FWRITEINT((int)NY(opts.cir));
				FWRITEARR(P(opts.cir), sizeof(real)*PN(opts.cir));
			}

			if(opts.limit){/*xmin,xmax,ymin,ymax */
				FWRITECMDARR(DRAW_LIMIT, opts.limit, sizeof(real)*4);
			}

			FWRITECMDSTR(DRAW_TITLE, title);
			FWRITECMDSTR(DRAW_XLABEL, xlabel);
			FWRITECMDSTR(DRAW_YLABEL, ylabel);

			FWRITEINT(DRAW_END);
end2:
			BUF_POSTPROC;
		}
		if(bufsize&&!ans){
			ans=send_buf(fig, fn, buf, bufsize, opts.always);
			free_default(buf);
		}
	}
	return ans;
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
/*
	Compute the graph limit. ox and oy are low left corner of the grid. When plotting, shown on the center of the the lower left "pixel"
	Seperate out to be easily modifiable.
*/
#define LIMIT_SET_X(limit, ox, offset, dx, nx) limit[0]=(ox)+(offset)*fabs(dx); limit[1]=limit[0]+fabs(dx)*(nx);
#define LIMIT_SET_Y(limit, ox, offset, dx, nx) limit[2]=(ox)+(offset)*fabs(dx); limit[3]=limit[2]+fabs(dx)*(nx);

/**
   Mapping the floating point numbers onto screen with scaling similar to matlab
   imagesc.  . see ddraw()
*/

/**
   like ddraw, acting on map object. see ddraw()

   @param[in] fig Category of the figure.
   @param[in] map OPD map to plot.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawmap(const char* fig, const map_t* map, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!map || !draw_current(fig, fn)) return 0;
	real limit[4];
	LIMIT_SET_X(limit, map->ox, 0.5, map->dx, map->nx);
	LIMIT_SET_Y(limit, map->oy, 0.5, map->dy, map->ny);
	draw(fig, (plot_opts){.image=(const dmat*)map, .zlim={-zlim,zlim}, .limit=limit}, title, xlabel, ylabel, "%s", fn);
	return 1;
}
/**
   Plot the loc on the screen. see ddraw()

   @param[in] fig Category of the figure.
   @param[in] loc Coordinate to plot. Its map is created if not yet.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawloc(const char* fig, loc_t* loc, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!loc || !draw_current(fig, fn)) return 0;
	loc_create_map(loc);
	int npad=loc->npad;
	long nx=loc->map->nx-npad*2;
	long ny=loc->map->ny-npad*2;
    dmat *opd0=dnew(nx,ny);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(opd0, ix, iy)=(P(loc->map, (ix+npad), (iy+npad))>0);
		}
	}
	real offset=npad+(isfinite(loc->ht)?-0.5:0);//-0.5 means coordinate is at center of pixel. 0 means at corner. saloc is at corner. 
	real limit[4];
	LIMIT_SET_X(limit, loc->map->ox, offset, loc->dx, nx);
	LIMIT_SET_Y(limit, loc->map->oy, offset, loc->dy, ny);
	draw(fig, (plot_opts){.image=opd0, .zlim={-zlim,zlim}, .limit=limit}, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}

/**
   Plot the opd using coordinate loc. see ddraw()

   @param[in] fig Category of the figure.
   @param[in] loc Coordinate of opd. Its map is created if not yet.
   @param[in] opd OPD defined on loc.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawopd(const char* fig, loc_t* loc, const dmat* opd, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){

	format2fn;
	if(!loc || !opd || !draw_current(fig, fn)) return 0;
	if(loc->nloc!=PN(opd)){
		warning("Invalid dimensions. loc has %ld, opd has %ldx%ld\n", loc->nloc, NX(opd), NY(opd));
		return 0;
	}
	dmat* opd0=dnew(0,0);
	loc_embed(opd0, loc, opd, 0, 0);
	real offset=loc->npad+(isfinite(loc->ht)?-0.5:0);//-0.5 means coordinate is at center of pixel. 0 means at corner. saloc is at corner. 
	real limit[4];
	if(loc->map){
		LIMIT_SET_X(limit, loc->map->ox, offset, loc->dx, opd0->nx);
		LIMIT_SET_Y(limit, loc->map->oy, offset, loc->dy, opd0->ny);
	}
	draw(fig, (plot_opts){.image=opd0, .zlim={-zlim,zlim}, .limit=limit}, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}
/**
   Sanitize gradients for plotting

   Set to NaN those where saa<0.01
   Remove tip/tilt if trs is set

   @param[in] gradin Gradients with all x gradients first.
   @param[in] saa Subaperture area. NULL to skip masking.
   @param[in] nsa Number of subapertures.
   @param[in] trs Whether to remove tip/tilt (only for 2 gradients per subaperture).
   @return Gradients with the same shape as gradin. A copy if modified, otherwise a reference. The caller frees it.
*/
static dmat* grad_prep(const dmat *gradin, const dmat *saa, long nsa, int trs){
	long ng=PN(gradin)/nsa;
	dmat *grad=0;
	if((trs && ng==2) || saa){
		grad=ddup(gradin);
	}else{
		grad=dref(gradin); 
	}
	reshape(grad, nsa, ng);
	if(saa){
		for(long isa=0; isa<nsa; isa++){
			if(P(saa, isa)<0.01){
				P(grad, isa, 0)=NAN;
				P(grad, isa, 1)=NAN;
			}
		}
	}
	if(trs&&ng==2){//remove tip/tilt
		real gxm=0;
		real gym=0;
		for(long isa=0; isa<nsa; isa++){
			if(!saa || P(saa, isa)>0.1){
				gxm+=P(grad, isa, 0);
				gym+=P(grad, isa, 1);
			}
		}
		gxm/=-nsa;
		gym/=-nsa;
		for(long isa=0; isa<nsa; isa++){
			if(!saa||P(saa, isa)>0.5){
				P(grad, isa, 0)+=gxm;
				P(grad, isa, 1)+=gym;
			}
		}
	}
	reshape(grad, NX(gradin), NY(gradin));
	return grad;
}
/**
   Plot gradients, either reconstructed into OPD using CuReD or as x and y gradient maps side by side.

   @param[in] fig Category of the figure.
   @param[in] saloc Subaperture locations.
   @param[in] saa Subaperture area used to mask gradients when grad2opd is set. May be NULL.
   @param[in] gradin Gradients with all x gradients first.
   @param[in] grad2opd Whether to reconstruct OPD with CuReD (requires more than 4 subapertures and 2 gradients each).
   @param[in] trs Whether to remove tip/tilt.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawgrad(const char* fig, loc_t* saloc, const dmat *saa, const dmat* gradin, int grad2opd, int trs, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!draw_current(fig, fn) || !saloc || !gradin) return 0;
	long nsa=saloc->nloc;
	long ng=PN(gradin)/nsa;
	if(nsa<=4||ng!=2) grad2opd=0;

	dmat *grad=grad_prep(gradin, grad2opd?saa:NULL, nsa, trs);
	dmat *phi=0;
	
	dmat *phix=dnew(0, 0);
	dmat *phiy=dnew(0, 0);
	dmat *gx=dnew_do(nsa, 1, P(grad), 0);
	dmat *gy=dnew_do(nsa, 1, P(grad)+nsa, 0);
	loc_embed(phix, saloc, gx, 0, 0);
	loc_embed(phiy, saloc, gy, 0, 0);
	dfree(gx);
	dfree(gy);
	if(grad2opd){
		cure(&phi, phix, phiy, saloc->dx);
	}else{
		phi=dcat(phix, phiy, 1);
	}
	dfree(phix);
	dfree(phiy);

	//This is different from loc_embed. It removes the padding.
	real offset=saloc->npad+(isfinite(saloc->ht)?-0.5:0);//-0.5 means coordinate is at center of pixel. 0 means at corner. saloc is at corner. 
	real limit[4];
	LIMIT_SET_X(limit, saloc->map->ox, offset, saloc->dx, phi->nx);
	LIMIT_SET_Y(limit, saloc->map->oy, offset, saloc->dy, phi->ny);
	draw(fig, (plot_opts){
		.image=phi, .zlim={-zlim,zlim}, .limit=limit}, title, xlabel, ylabel, "%s", fn);
	
	dfree(grad);
	dfree(phi);
	return 1;
}
/**
   Plot opd with coordinate loc where amp is above half of its maximum. see ddraw()

   @param[in] fig Category of the figure.
   @param[in] loc Coordinate of opd. Its map is created if not yet.
   @param[in] opd OPD defined on loc.
   @param[in] amp Amplitude defined on loc.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawopdamp(const char* fig, loc_t* loc, const dmat* opd, const dmat* amp, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!loc || !opd || !amp || !draw_current(fig, fn)) return 0;
	(void)fig;
	if(loc->nloc!=NX(amp)||loc->nloc!=PN(opd)){
		warning("Invalid dimensions. loc has %ld, opd has %ldx%ld, amp has %ldx%ld.\n",
			loc->nloc, NX(opd), NY(opd), NX(amp), NY(amp));
		return 0;
	}
	loc_create_map(loc);

	int npad=loc->npad;
	long nx=loc->map->nx-npad*2;
	long ny=loc->map->ny-npad*2;
	real ampthres;
	dvecmaxmin(P(amp), loc->nloc, &ampthres, 0);
	ampthres*=0.5;
    dmat* opd0=dnew(nx,ny);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			long ii=P(loc->map, (ix+npad), (iy+npad))-1;
            P(opd0, ix, iy)=(ii>-1&&P(amp, ii)>ampthres)?P(opd, ii):NAN;
		}
	}
	real offset=npad+(isfinite(loc->ht)?-0.5:0);//-0.5 means coordinate is at center of pixel. 0 means at corner. saloc is at corner. 
	real limit[4];
	LIMIT_SET_X(limit, loc->map->ox, offset, loc->dx, nx);
	LIMIT_SET_Y(limit, loc->map->oy, offset, loc->dy, ny);
	draw(fig, (plot_opts){
		.image=opd0, .zlim={-zlim,zlim}, .limit=limit}, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}
/**
   Concatenate and plot subaperture images.

   @param[in] fig Category of the figure.
   @param[in] saloc Subaperture locations.
   @param[in] ints Subaperture images: one per subaperture, or a single (TT or PWFS) or 4 (TTF) images.
   @param[in] zlim Limit of the color scale is [-zlim, zlim]. 0 for automatic.
   @param[in] title Title of the plot.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the plot name, followed by its arguments.
   @return 1 if the plot is sent to draw(), 0 if skipped (not current page or invalid input).
*/
int drawints(const char* fig, const loc_t* saloc, const dcell* ints, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!saloc || !ints || !draw_current(fig, fn)) return 0;
	dmat* ints2=0;
	if(NX(ints)==1){//TT or PWFS
		if(P(ints, 0)->nx==P(ints, 0)->ny){//TT
			ints2=dref(P(ints, 0));
		} else{//PWFS
			dcell* ints3=loc_embed2(saloc, P(ints, 0));
			if(NX(ints3)*NY(ints3)==4){
				reshape(ints3, 2, 2);
			}
			ints2=dcell2m(ints3);
			dcellfree(ints3);
		}
	} else if(NX(ints)==4){//TTF
		dcell* ints3=dcellref(ints);
		reshape(ints3, 2, 2);
		ints2=dcell2m(ints3);
		dcellfree(ints3);
	} else{
		dcell* ints3=0;
		loc_embed_cell(&ints3, saloc, ints);
		ints2=dcell2m(ints3);
		dcellfree(ints3);
	}
	draw(fig, (plot_opts){
		.image=ints2, .zlim={-zlim,zlim}}, title, xlabel, ylabel, "%s", fn);
	dfree(ints2);
	return 1;
}
