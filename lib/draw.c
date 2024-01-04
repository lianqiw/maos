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
PNEW(lock);
#define MAXDRAW 1024
#define TEST_UDP 0 //test UDP implementation. Doesn't seem to help performance. Keep at 0.

int draw_id=0;      //Client identification. Same draw_id  reuses drawdaemon.
int draw_direct=0;  //Directly launch drawdaemon without forking a draw_helper
int draw_disabled=0; //if set, draw will be disabled
static int sock_helper=-1;//socket of draw_helper
int listening=0;    //If set, listen_drawdaemon is listening to replies from drawdemon
int draw_single=0;  //if ==1, only draw active frame and skip if line is busy. if ==-1, disable draw_single altogether. otherwise draw all frames.

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
   draw_helper() will then create another socket_pair, pass one end of it to
   draw(), and the other end to drawdaemon() by fork()+exec().

   3 If draw_direct=1, get_drawdaemon() will launch drawdaemon directly.

*/

/* List of list*/
typedef struct list_t{
	char* key;
	struct list_t* next;
	struct list_t* child;//child list
}list_t;

typedef struct sockinfo_t{
	int fd;//socket;
#if TEST_UDP
	udp_t udp;
#endif
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
   Listen to drawdaemon for update of fig, fn. The hold values are stored in figfn.
*/
static void* listen_drawdaemon(sockinfo_t* sock_data){
	listening=1;
	int sock_draw=sock_data->fd;
	char** figfn=sock_data->figfn;
	//info("draw is listening to drawdaemon at %d\n", sock_draw);
#if TEST_UDP	
	if(sock_data->udp.sock<=0){
		sock_data->udp.sock=bind_socket(SOCK_DGRAM, 0, 0);
		int client_port=socket_port(sock_data->udp.sock);
		int cmd[4]={DRAW_ENTRY, sizeof(int)*2, DRAW_UDPPORT, client_port};
		static int nretry=0;
		dbg("client_port=%d\n", client_port);
retry:
		if(stwrite(sock_draw, cmd, sizeof(cmd))){
			warning("write to drawdaemon failed\n");
			return NULL;
		}
		//server replies in UDP packet.
		socklen_t slen=sizeof(sock_data->udp.peer_addr);	
		int cmd2[64];
		socket_recv_timeout(sock_data->udp.sock, 2);
		int ncmd=recvfrom(sock_data->udp.sock, cmd2, sizeof(cmd2), 0, 
						 (struct sockaddr*)&sock_data->udp.peer_addr, &slen);
		if(ncmd>0){
			dbg("received udp reply from %s:%d (%d) bytes. UDP payload is %d. Header size is %d\n", 
				addr2name(sock_data->udp.peer_addr.sin_addr.s_addr), 
						  sock_data->udp.peer_addr.sin_port, cmd2[3], cmd2[4], cmd2[5]);
			sock_data->udp.version=cmd2[2];
			sock_data->udp.payload=cmd2[4];
			sock_data->udp.header=cmd2[5];
			//set default destination for UDP.
			if(connect(sock_data->udp.sock, (struct sockaddr*)&sock_data->udp.peer_addr, slen)){
				dbg("connection fails to establish, errno=%d, close socket.\n", errno);
				close(sock_data->udp.sock);
				sock_data->udp.sock=-1;
			}else{
				dbg("connection is established.\n");
			}
		}else{
			nretry++;
			if(nretry<5){
				warning("no data is received, retry\n");
				goto retry;
			}else{
				warning("no data is received, give up\n");
			}
		}
	}
#endif	
	int cmd;
	int nlen=0;

	if(!socket_recv_timeout(sock_draw, 0)){
        //make sure it blocks when no data is readable
        dbg("started listening at %d\n", sock_draw);
    }else{
        dbg("set timeout=0 failed for %d\n", sock_draw);
    }
	while(sock_data->fd!=-1&&!streadint(sock_draw, &cmd)){
		if(cmd==DRAW_ENTRY){//every message in new format start with DRAW_ENTRY.
			streadint(sock_draw, &nlen);
			streadint(sock_draw, &cmd);
		}
		switch(cmd){
		case DRAW_FIGFN:
		{
			char* fig=0, * fn=0;
			streadstr(sock_draw, &fig);
			streadstr(sock_draw, &fn);
			if(fig&&fn){
				if(figfn[0]&&figfn[1]&&(strcmp(figfn[0], fig)||strcmp(figfn[1], fn))){
					if(sock_data->draw_single==1){
						dbg2("draw %d switch to fig=%s, fn=%s\n", sock_draw, fig, fn);
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
			dbg2("draw %d paused\n", sock_draw);
			break;
		case DRAW_RESUME:
			sock_data->pause=0;
			dbg2("draw %d resumed\n", sock_draw);
			break;
		case DRAW_SINGLE:
			sock_data->draw_single=!sock_data->draw_single;
			dbg2("draw %d draw_single=%d\n", sock_draw, sock_data->draw_single);
			break;
		default:
			dbg("Unknown cmd: %d with size %d from socket %d\n", cmd, nlen, sock_draw);
			if(nlen){
				void* p=malloc(nlen);
				stread(sock_draw, p, nlen);
				free(p);
			}
			break;
		}
	}
	dbg_time("stopped lisening to drawdaemon at %d, errno=%d, %s\n", sock_draw, errno, strerror(errno));
	listening=0;
	return NULL;
}

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
/**Add fd to list of drawing socks*/
int draw_add(int sock){
	//Check that the drawdaemon is live.
	if(sock!=-1){
		int cmd[2]={DRAW_PID, getpid()};
		if(stwrite(sock, cmd, sizeof(cmd))){
			dbg("failed to write: sock=%d\n", sock);
			close(sock); sock=-1;
		}
	}
	if(sock==-1) return -1;
	sockinfo_t *p=mycalloc(1, sockinfo_t);
	p->fd=sock;
	p->next=sock_draws;
	sock_draws=p;
	if(draw_single!=-1){
		p->draw_single=1;
	    pthread_create(&p->thread, NULL, (thread_fun)listen_drawdaemon, p);
	}
	draw_disabled=0;
	return 0;
}
///fd==-1 will remove all clients
///reuse==0 is called when write failed, do not try to write again
static void draw_remove(int fd, int reuse){
	for(sockinfo_t **curr=&sock_draws; *curr;){
		sockinfo_t *p=*curr;
		if(p->fd==fd || fd==-1){
			*curr=p->next;
			if(reuse){
				stwriteint(p->fd, DRAW_FINAL);
				dbg("send %d back to scheduler for reuse\n", p->fd);
				scheduler_socket(1, &p->fd, draw_id?draw_id:1);
			}
			if(p->thread){
				void *ans;
				if(pthread_cancel(p->thread) || pthread_join(p->thread, &ans)){
					dbg("Unable to cancel or join thread\n");
				}
			}
			if(fd!=-1) close(fd);
#if TEST_UDP			
			if(p->udp.sock>0) {
				close(p->udp.sock);
			}
#endif
			list_destroy(&p->list);
			free(p->figfn[0]);
			free(p->figfn[1]);
			free(p);
		}else{
			curr=&p->next;
		}
	}
	if(!sock_draws){
		draw_disabled=1;//do not try to restart drawdaemon
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
	} else{
		perror("socketpair");
		warning("socket pair failed, cannot launch drawdaemon\n");
		draw_disabled=1;
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
		warning("socketpair failed, disable drawing.\n");
		draw_disabled=1;
	} else{
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
	} else{/*Child */
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
	if(sock_draws){//drawdaemon already connected
		return 0;
	}
	if(draw_disabled){
		return -1;
	}
#if __APPLE__
	int display=1;
#else
	char* display=getenv("DISPLAY");
	if(display&&!strlen(display)){//display is not set, we ask scheduler to open a drawdaemon.
		display=0;
	}
#endif
	if(draw_id<=0){
		draw_id=getsid(0);
		if(draw_id<=0){
			draw_id=1;
		}
	}
	int sock=-1;
	//First try reusing existing idle drawdaemon with the same id
	while(!scheduler_socket(-1, &sock, draw_id)){
		//test whether received drawdaemon is still running
		if(stwriteint(sock, DRAW_INIT)){
			dbg("received socket=%d is already closed.\n", sock);
			close(sock);
			sock=-1;
		}else{
			dbg("received socket=%d is ok.\n", sock);
			break;
		}
	}
	if(sock==-1){
		if(display){
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
					warning("Unable to talk to the helper to launch drawdaemon.\n");
				}else{
					dbg("launch using sock helper: sock=%d\n", sock);
				}
			}
		} else{//no display is available. use scheduler to launch drawdaemon
			scheduler_socket(-1, &sock, 0);
			dbg("launch using scheduler: sock=%d\n", sock);
		}
	}

	if(sock!=-1 && !draw_add(sock) && !socket_send_timeout(sock, 60)){
		//prevent hang. too small timeout will prevent large data from passing through.
		return 0;
	}else{
		draw_disabled=1;
		warning("Unable to open drawdaemon. Disable drawing.\n");
		return 1;
	}
}
/**
   Tell drawdaemon that this client will no long use the socket. Send the socket to scheduler for future reuse.
*/
void draw_final(int reuse){
	//called from other threads, need to lock
	LOCK(lock);
	draw_remove(-1, reuse);
	UNLOCK(lock);
}

/*
   Check whether we need to draw current page. True if
   1) not yet exist in the list. it is inserted into the list if add is valid and set as current page.
   2) is current page.
   3) when draw_single is not set.
   4) pause must not be set set.

   return
   1: always draw
   2: draw if lock success
   3: create an empty page.
*/
static int check_figfn(sockinfo_t *ps, const char* fig, const char* fn){
	if(draw_disabled||ps->fd==-1||ps->pause) return 0;
	if(!(draw_single==1 && ps->draw_single==1)) return 1;
	list_t* child=0;
	list_search(&ps->list, &child, fig, 1);
	int found=0;
	if(child){
		if(fn){
			found=list_search(&child->child, NULL, fn, 1);
		} else{//don't check fn
			found=1;
		}
	}
	int ans=0;//default is false
	if(!found){//page is not found
		//ans=3; //draw without test lock
		plot_empty(ps->fd, fig, fn);
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
	return ans;
}

/**
   Check whether what we are drawing is current page of any drawdaemon.
*/
int draw_current(const char* fig, const char* fn){
    int current=0;
	if(draw_disabled){
        current=0;
    }else if(draw_single!=1){
        current=1;
    }else{
		for(sockinfo_t *ps=sock_draws; ps; ps=ps->next){
            if((current=check_figfn(ps, fig, fn))){
                break;
            }
        }
        if(current==2){//check whether draw is busy
			static int nskip=0;
			static int nplot=0;
            if((TRYLOCK(lock))){//lock failed
                current=0;
				nskip++;
                if(nplot) dbg2("Skip after plot %d times\n", nplot);
				nplot=0;			
            } else{
				nplot++;
				if(nskip) dbg2("Resume after skip %d times\n", nskip);
				nskip=0;
                UNLOCK(lock);
            }
        }
    }
	return current;
}
static inline int fwriteint(FILE* fbuf, int A){
	if(fwrite(&A, sizeof(int), 1, fbuf)<1) return -1;
	return 0;
}
//write str length and str itself (with tailing NULL)
static inline int fwritestr(FILE* fbuf, const char* str){
	if(!str) return 0;
	uint32_t nlen=strlen(str)+1;
	if(fwriteint(fbuf, nlen)||fwrite(str, 1, nlen, fbuf)<nlen) return -1;
	return 0;
}
static int iframe=0; //frame counter
#define CATCH(A) if(A) {ans=1; warning("fwrite failed\n"); goto end2;}
#define FWRITEARR(p,len) CATCH(fwrite(p,1,len,fbuf)<len)
#define FWRITESTR(str) CATCH(fwritestr(fbuf, str))
#define FWRITEINT(A) CATCH(fwriteint(fbuf,A))
#define FWRITECMD(cmd, nlen) CATCH(fwriteint(fbuf, DRAW_ENTRY) || fwriteint(fbuf, (nlen)+sizeof(int)) || fwriteint(fbuf, cmd))
#define FWRITECMDSTR(cmd,str) if(str){FWRITECMD(cmd, strlen(str)+1); FWRITESTR(str);}
#define FWRITECMDARR(cmd,p,len) if((len)>0){FWRITECMD(cmd, len); FWRITEARR(p,len);}
//Plot an empty page
int plot_empty(int sock_draw, const char *fig, const char *fn){
	char *buf=0;
	int ans=0;
	size_t bufsize=0;
	FILE *fbuf=open_memstream(&buf,&bufsize);
	int zeros[4]={0,0,0,0};
	FWRITECMDARR(DRAW_FRAME,zeros,4*sizeof(int));
	FWRITECMD(DRAW_START,0);
	FWRITECMD(DRAW_FLOAT,sizeof(int));FWRITEINT(sizeof(real));
	FWRITECMDSTR(DRAW_FIG,fig);
	FWRITECMDSTR(DRAW_NAME,fn);
	FWRITECMD(DRAW_END,0);
end2:
	fclose(fbuf);
	if(stwrite(sock_draw,buf,bufsize)){
		info("write to %d failed: %s\n",sock_draw,strerror(errno));
		ans=-1;
		draw_remove(sock_draw,0);
	}
	free_default(buf);
	return ans;
}
int send_buf(const char *fig, const char *fn, char *buf, size_t bufsize){
	int ans=0;
	for(sockinfo_t *ps=sock_draws, *ps_next=NULL; ps; ps=ps_next){
		ps_next=ps->next;//draw_remove change ps node, so save info here;
		/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
		int sock_draw=ps->fd;
		int needed_i=check_figfn(ps, fig, fn);
		if(!needed_i){
			continue;
		}
		if(needed_i==2){
			if(TRYLOCK(lock)){//line busy
				continue;
			}
		} else{
			LOCK(lock);
		}
		if(use_udp){
			error("To be implemented\n");
			//use sendmmsg with GSO is fastest.
		} else{
			//TIC;tic;
			if(stwrite(sock_draw, buf, bufsize)){
				info("write to %d failed: %s\n", sock_draw, strerror(errno));
				ans=-1;
				draw_remove(sock_draw, 0);
			} else{
				//toc2("write %s:%s %lu MiB", fig, fn, bufsize>>20);
			}
		}
		
#if TEST_UDP				
		if(ps->udp.sock>0){
			//test UDP data
			int counter=myclocki();
			udp_send(&ps->udp, buf, bufsize, counter);
			udp_send(&ps->udp, buf, bufsize, counter);
			udp_send(&ps->udp, buf, bufsize, counter);
		}
#endif	
		UNLOCK(lock);
	}
	return ans;
}
/**
   Plot the coordinates ptsx, ptsy using style, and optionally plot ncir circles.
*/
int plot_points(const char* fig,    /**<Category of the figure*/
	plot_opts opts,     /**<Additional options*/
	const char* title,  /**<title of the plot*/
	const char* xlabel, /**<x axis label*/
	const char* ylabel, /**<y axis label*/
	const char* format, /**<subcategory of the plot.*/
	...){
	if(draw_disabled) return 0;
	format2fn;
	
	iframe++;
	int ans=0;
	if(!get_drawdaemon()){
		int needed=0;
		for(sockinfo_t *ps=sock_draws; ps; ps=ps->next){
			/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
			int needed_i=check_figfn(ps, fig, fn);
			//dbg("%s %s: %d\n", fig, fn, needed_i);
			if(needed_i){
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
			int zeros[4]={0,0,0,0};	
			FWRITECMDARR(DRAW_FRAME, zeros, 4*sizeof(int));
			FWRITECMD(DRAW_START, 0);
			FWRITECMD(DRAW_FLOAT, sizeof(int));FWRITEINT(sizeof(real));
			FWRITECMDSTR(DRAW_FIG, fig);
			FWRITECMDSTR(DRAW_NAME, fn);
			if(opts.loc){/*there are points to plot. */
				for(int ig=0; ig<opts.ngroup; ig++){
					int nlen=opts.loc[ig]->nloc;
					if(opts.maxlen && opts.maxlen<nlen) nlen=opts.maxlen;
					FWRITECMD(DRAW_POINTS, 3*sizeof(int)+sizeof(real)*nlen*2);
					FWRITEINT(nlen);
					FWRITEINT(2);
					FWRITEINT(1);
					FWRITEARR(opts.loc[ig]->locx, sizeof(real)*nlen);
					FWRITEARR(opts.loc[ig]->locy, sizeof(real)*nlen);
				}
				if(opts.dc){
					warning("both loc and dc are specified, ignore dc.\n");
				}
			} else if(opts.dc){
				if(opts.ngroup>PN(opts.dc)||opts.ngroup==0){
					opts.ngroup=PN(opts.dc);
				}
				for(int ig=0; ig<opts.ngroup; ig++){
					dmat* p=P(opts.dc, ig);
					int nlen=NX(p);
					if(opts.maxlen&&opts.maxlen<nlen) nlen=opts.maxlen;
					FWRITECMD(DRAW_POINTS, 3*sizeof(int)+nlen*sizeof(real));
					FWRITEINT(nlen);//number of points
					FWRITEINT(NY(p));//number of numbers per point. 1 or 2.
					FWRITEINT(0);//square plot or not
					if(p){
						FWRITEARR(P(p), NY(p)*nlen*sizeof(real));
					}
				}
			} else{
				error("Invalid Usage\n");
			}
			if(opts.style){
				FWRITECMD(DRAW_STYLE, sizeof(int)+sizeof(uint32_t)*opts.ngroup);
				FWRITEINT(opts.ngroup);
				FWRITEARR(opts.style, sizeof(uint32_t)*opts.ngroup);
			}
			if(opts.cir){
				if(NX(opts.cir)!=4){
					error("Cir should have 4 rows\n");
				}
				FWRITECMD(DRAW_CIRCLE, sizeof(int)+sizeof(real)*PN(opts.cir));
				FWRITEINT(NY(opts.cir));
				FWRITEARR(P(opts.cir), sizeof(real)*PN(opts.cir));
			}
			if(opts.limit){/*xmin,xmax,ymin,ymax */
				FWRITECMDARR(DRAW_LIMIT, opts.limit, sizeof(real)*4);
			}
			if(opts.xylog){
				FWRITECMDARR(DRAW_XYLOG, opts.xylog, sizeof(char)*2);
			}
			/*if(format){
				FWRITECMDSTR(DRAW_NAME, fn);
			}*/
			if(opts.legend){
				int nlen=0;
				for(int ig=0; ig<opts.ngroup; ig++){
					if(opts.legend[ig]){
						nlen+=strlen(opts.legend[ig])+1;
					}
				}
				FWRITECMD(DRAW_LEGEND, nlen);
				for(int ig=0; ig<opts.ngroup; ig++){
					FWRITESTR(opts.legend[ig]);
				}
			}

			FWRITECMDSTR(DRAW_TITLE, title);
			FWRITECMDSTR(DRAW_XLABEL, xlabel);
			FWRITECMDSTR(DRAW_YLABEL, ylabel);
		
			FWRITECMD(DRAW_END, 0);
end2:
			fclose(fbuf);
			if(ans){
				warning("write failed:%d\n", ans);
				free_default(buf); bufsize=0; buf=0;
			} else{
				int* bufp=(int*)(buf+3*sizeof(int));
				bufp[0]=(int)bufsize;//frame size
				bufp[1]=iframe;//frame number
				bufp[2]=(int)bufsize;//sub-frame size
				bufp[3]=0;//sub-frame number
			}
		}
		if(bufsize&&!ans){
			ans=send_buf(fig, fn, buf, bufsize);
			free_default(buf);
		}
	}

	return ans;
}

/**
   Draw an image.
*/
//Data structure for imagesc
#define MAXNX 1000 //downsample
typedef float dtype;
typedef struct imagesc_t{
	char* fig; /**<Category of the figure*/
	long nx;   /**<the image is of size nx*ny*/
	long ny;   /**<the image is of size nx*ny*/
	long bsize;/**<bytes of each element*/
	dtype* limit; /**<x min; xmax; ymin and ymax*/
	dtype* zlim;  /**<Clip the data to min; max*/
	dtype* p; 	  /**<The image*/
	char* title;  /**<title of the plot*/
	char* xlabel; /**<x axis label*/
	char* ylabel; /**<y axis label*/
	char* fn;
}imagesc_t;
static void* imagesc_do(imagesc_t* data){
	if(draw_disabled) return NULL;
	iframe++;//frame counter
	
	if(!get_drawdaemon()){
		char* fig=data->fig;
		long nx=NX(data);
		long ny=NY(data);
		dtype* limit=data->limit;
		dtype* zlim=data->zlim;
		dtype* p=P(data);
		char* title=data->title;
		char* xlabel=data->xlabel;
		char* ylabel=data->ylabel;
		char* fn=data->fn;
		int needed=0;
		for(sockinfo_t *ps=sock_draws; ps; ps=ps->next){
			/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
			int needed_i=check_figfn(ps, fig, fn);
			//dbg("%s %s: %d\n", fig, fn, needed_i);
			if(needed_i){
				if(!needed||needed==3){
					needed=needed_i;
				}
			}
		}
		int ans=0;
		char* buf=0;
		size_t bufsize=0;
		if(needed){
			FILE* fbuf=open_memstream(&buf, &bufsize);
			int zeros[4]={0,0,0,0};
			FWRITECMDARR(DRAW_FRAME, zeros, 4*sizeof(int));
			FWRITECMD(DRAW_START, 0);
            FWRITECMD(DRAW_FLOAT, sizeof(int));FWRITEINT(sizeof(dtype));
			FWRITECMDSTR(DRAW_FIG, fig);
			FWRITECMDSTR(DRAW_NAME, fn);
			
			size_t nlen=2*sizeof(int32_t)+sizeof(dtype)*nx*ny;
			FWRITECMD(DRAW_DATA, nlen);
			int32_t header[2];
			header[0]=nx;
			header[1]=ny;
			FWRITEARR(header, sizeof(int32_t)*2);
			FWRITEARR(p, sizeof(dtype)*nx*ny);
			if(zlim){
				FWRITECMDARR(DRAW_ZLIM, zlim, sizeof(dtype)*2);
			}
			if(limit){/*xmin,xmax,ymin,ymax */
				FWRITECMDARR(DRAW_LIMIT, limit, sizeof(dtype)*4);
			}

			FWRITECMDSTR(DRAW_TITLE, title);
			FWRITECMDSTR(DRAW_XLABEL, xlabel);
			FWRITECMDSTR(DRAW_YLABEL, ylabel);
			
			FWRITECMD(DRAW_END, 0);
end2:
			fclose(fbuf);
			if(ans){
				warning("write failed:%d\n", ans);
				free_default(buf); bufsize=0; buf=0;
			} else{
				int* bufp=(int*)(buf+3*sizeof(int));
				bufp[0]=(int)bufsize;//frame number
				bufp[1]=iframe;//frame number
				bufp[2]=(int)bufsize;//sub-frame number
				bufp[3]=0;//sub-frame number
				//dbg("draw_image: buf has size %ld. %d %d %d %d\n", bufsize,
				//	bufp[0], bufp[1], bufp[2], bufp[3]);
			}
		}

		if(bufsize&&!ans){
			ans=send_buf(fig, fn, buf, bufsize);
			free_default(buf);
		}
	}
	free(data->fig);
	free(data->limit);
	free(data->zlim);
	free(P(data));
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
int ddraw(const char *fig, /**<Category of the figure*/
            const dmat *p,   /**<The image*/
            const real *limit,  /**<x min, xmax, ymin and ymax*/
            const real *zlim,   /**< min,max, the data*/
            const char *title,  /**<title of the plot*/
            const char *xlabel, /**<x axis label*/
            const char *ylabel, /**<y axis label*/
            const char *format, /**<subcategory of the plot.*/
            ...){
    format2fn;
    if(!p||!draw_current(fig, fn)){
        return 0;
    }

    // We copy all the data and put the imagesc job into a task
    // The task will free the data after it finishes.
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
        int xstep=(NX(p)+MAXNX-1)/MAXNX;
        int ystep=(NY(p)+MAXNX-1)/MAXNX;
        int nx2=(NX(p))/xstep;
        int ny2=(NY(p))/ystep;
        P(data)=malloc(nx2*ny2*sizeof(dtype));
        for(int iy=0; iy<ny2; iy++){
            dtype *p2=P(data)+iy*nx2;
            const real *p1=P(p)+iy*ystep*NX(p);
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
	//dbg("draw_single=%d\n", draw_single);
    if(draw_single==1){
        thread_new((thread_fun)imagesc_do, data);
    } else{
        imagesc_do(data);
    }
    return 1;
}

/**
   Draw the OPD of abs and phase of complex p defined on nx*ny grid. see ddraw()
   type: 0: abs only. 1: abs and phase. 2: real and imaginary
*/
int cdraw(const char *fig, const cmat *p, const real *limit, const real *zlim,
    int type, const char *title, const char *xlabel, const char *ylabel,
    const char *format, ...){
    format2fn;
    if(!draw_current(fig, fn)) return 0;
    int isreal=1;
    long nx=NX(p);
    long ny=NY(p);
    dmat *pr=dnew(nx, ny);
    dmat *pi=dnew(nx, ny);
    for(int i=0; i<nx*ny; i++){
        comp ppi=P(p, i);
        if(type==2){
            P(pr, i)=creal(ppi);
            P(pi, i)=cimag(ppi);
        } else{
            P(pr, i)=cabs(ppi);
            if(type==1&&P(pr, i)>1.e-10){
                P(pi, i)=atan2(cimag(ppi), creal(ppi));
            }
        }
        if(isreal&&fabs(P(pi, i))>1.e-10) isreal=0;
    }
    ddraw(fig, pr, limit, zlim, title, xlabel, ylabel, "%s abs", fn);
    dfree(pr);
    if(!isreal){
        ddraw(fig, pi, limit, zlim, title, xlabel, ylabel, "%s phi", fn);
    }
    dfree(pi);
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
*/
int drawmap(const char* fig, const map_t* map, real* zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!map || !draw_current(fig, fn)) return 0;
	real limit[4];
	LIMIT_SET_X(limit, map->ox, 0.5, map->dx, map->nx);
	LIMIT_SET_Y(limit, map->oy, 0.5, map->dx, map->ny);
	/*limit[0]=map->ox-map->dx/2;
	limit[1]=map->ox+(map->nx-0.5)*map->dx;
	limit[2]=map->oy-map->dx/2;
	limit[3]=map->oy+(map->ny-0.5)*map->dx;*/
	ddraw(fig, (const dmat*)map, limit, zlim, title, xlabel, ylabel, "%s", fn);
	return 1;
}
/**
   Plot the loc on the screen. see ddraw()
*/
int drawloc(const char* fig, loc_t* loc, real* zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!loc || !draw_current(fig, fn)) return 0;
	loc_create_map(loc);
	int npad=loc->npad;
	int nx=loc->map->nx-npad*2;
	int ny=loc->map->ny-npad*2;
    dmat *opd0=dnew(nx,ny);
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(opd0, ix, iy)=(P(loc->map, (ix+npad), (iy+npad))>0);
		}
	}
	real limit[4];
	LIMIT_SET_X(limit, loc->map->ox, npad-0.5, loc->dx, nx);
	LIMIT_SET_Y(limit, loc->map->oy, npad-0.5, loc->dx, ny);
	/*limit[0]=loc->map->ox+fabs(loc->dx)*(npad);
	limit[1]=limit[0]+loc->dx*(nx);
	limit[2]=loc->map->oy+fabs(loc->dy)*(npad);
	limit[3]=limit[2]+loc->dy*(ny);*/
	ddraw(fig, opd0, limit, zlim, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}

/**
   Plot the opd using coordinate loc. see ddraw()
*/
int drawopd(const char* fig, loc_t* loc, const dmat* opd, real* zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){

	format2fn;
	if(!loc || !opd || !draw_current(fig, fn)) return 0;
	if(loc->nloc!=PN(opd)){
		warning("Invalid dimensions. loc has %ld, opd has %ldx%ld\n", loc->nloc, NX(opd), NY(opd));
		return 0;
	}
	if(!loc->map){
		loc_create_map(loc);
	}
	//This is different from loc_embed. It removes the padding.
	int npad=loc->npad;
	int nx=loc->map->nx-npad*2;
	int ny=loc->map->ny-npad*2;
	dmat* opd0=dnew(nx, ny);
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			long ii=P(loc->map, (ix+npad), (iy+npad));
			P(opd0, ix, iy)=ii>0?P(opd, ii-1):NAN;
		}
	}
	real limit[4];
	LIMIT_SET_X(limit, loc->map->ox, npad-0.5, loc->dx, nx);
	LIMIT_SET_Y(limit, loc->map->oy, npad-0.5, loc->dy, ny);
/*	limit[0]=loc->map->ox+fabs(loc->dx)*(npad);
	limit[1]=limit[0]+fabs(loc->dx)*(nx);
	limit[2]=loc->map->oy+fabs(loc->dy)*(npad);
	limit[3]=limit[2]+fabs(loc->dy)*(ny);
	*/
	ddraw(fig, opd0, limit, zlim, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}
/**
   Plot gradients using CuReD
*/
int drawgrad(const char* fig, loc_t* saloc, const dmat* gradin, int grad2opd, int trs, real* zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(!saloc || !gradin) return 0;
    long nsa=saloc->nloc;
	if(nsa<=4||NY(gradin)!=2) grad2opd=0;
    //check current plotting target
    
	char fnx[100];
	char fny[100];
	if(NY(gradin)>1){
		snprintf(fnx, sizeof(fnx), "%sx", fn);
		snprintf(fny, sizeof(fny), "%sy", fn);
	}else{
		snprintf(fnx, sizeof(fnx), "%s", fn);
		snprintf(fny, sizeof(fny), "%s", fn);
	}
	if(!(grad2opd&&draw_current(fig, fn)) && !draw_current(fig, fnx) && !draw_current(fig, fny)) return 0;
	
	dmat* grad=0;
	if(trs && NY(grad)==2){//remove tip/tilt
		grad=ddup(gradin);
		reshape(grad, nsa, 2);
		real gxm=0;
		real gym=0;
		for(int isa=0; isa<nsa; isa++){
			gxm+=P(grad, isa, 0);
			gym+=P(grad, isa, 1);
		}
		gxm/=-nsa;
		gym/=-nsa;
		for(int isa=0; isa<nsa; isa++){
			P(grad, isa, 0)+=gxm;
			P(grad, isa, 1)+=gym;
		}
		reshape(grad, NX(gradin), NY(gradin));
	}else{
		grad=dref(gradin);
	}
	if(grad2opd&&draw_current(fig, fn)){
		//This is different from loc_embed. It removes the padding.
		dmat* phi=0;
		cure_loc(&phi, grad, saloc);
		real limit[4];
		int npad=saloc->npad;
		LIMIT_SET_X(limit, saloc->map->ox, npad-0.5, saloc->dx, phi->nx);
		LIMIT_SET_Y(limit, saloc->map->oy, npad-0.5, saloc->dy, phi->ny);
		ddraw(fig, phi, limit, zlim, title, xlabel, ylabel, "%s", fn);
		dfree(phi);
	}
	if(draw_current(fig, fnx)){
		dmat* gx=dnew_do(nsa, 1, P(grad), 0);
		drawopd(fig, saloc, gx, zlim, title, xlabel, ylabel, "%s", fnx);
		dfree(gx);
	}
	if(NY(grad)>1 && draw_current(fig, fny)){
		dmat *gy=dnew_do(nsa, 1, P(grad)+nsa, 0);
		drawopd(fig, saloc, gy, zlim, title, xlabel, ylabel, "%s", fny);
		dfree(gy);
	}
	dfree(grad);
	return 1;
}
/**
   Plot opd with coordinate loc where amp is above threshold. see ddraw()
*/
int drawopdamp(const char* fig, loc_t* loc, const dmat* opd, const dmat* amp, real* zlim,
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
	int nx=loc->map->nx-npad*2;
	int ny=loc->map->ny-npad*2;
	real ampthres;
	dmaxmin(P(amp), loc->nloc, &ampthres, 0);
	ampthres*=0.5;
    dmat* opd0=dnew(nx,ny);
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			long ii=P(loc->map, (ix+npad), (iy+npad))-1;
            P(opd0, ix, iy)=(ii>-1&&P(amp, ii)>ampthres)?P(opd, ii):NAN;
		}
	}
	real limit[4];
	LIMIT_SET_X(limit, loc->map->ox, npad-0.5, loc->dx, nx);
	LIMIT_SET_Y(limit, loc->map->oy, npad-0.5, loc->dy, ny);
	/*limit[0]=loc->map->ox+loc->dx*(npad-1);
	limit[1]=limit[0]+loc->dx*nx;
	limit[2]=loc->map->oy+loc->dx*(npad-1);
	limit[3]=limit[2]+loc->dx*ny;*/
	ddraw(fig, opd0, limit, zlim, title, xlabel, ylabel, "%s", fn);
	dfree(opd0);
	return 1;
}
/**
   Concatenate and plot subaperture images.
 */
int drawints(const char* fig, const loc_t* saloc, const dcell* ints, real* zlim,
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
	ddraw("Ints", ints2, NULL, zlim, title, xlabel, ylabel, "%s", fn);
	dfree(ints2);
	return 1;
}