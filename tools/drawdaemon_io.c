/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <errno.h>
#include "drawdaemon.h"
/*
  Routines in this file handles I/O.

  Todo: fread does not block when there are no more data available, and simply
  return EOF. Consider changing to read, which blocks when no data available.
 */
 //minimum internet maximum MTU is 576. IP header is 20-60. UDP header os 8. 508 is safest. We round up to 512.
#define UDP_PAYLOAD 512 //maximum size of UDP payload in bytes
#define UDP_HEADER 12 //size of UDP sub frame header in bytes
int ndrawdata=0;
int count=0;
int byte_float=sizeof(float);//default to float until client changes it.
udp_t udp_client={0};
int client_port=-1;//client udp port
in_addr_t client_addr;
int udp_sock=-1;//server udp socket
float io_time1=0;//time the latest drawdata is receievd
float io_time2=0;//time the previous drawdata is receveid
float io_timeclear=0;//plots
int io_heartbeat=0;
int session=1;//session counter
int noellipsis=1; 	/*do not allow legend ellipsis.*/
//This file does not link to math folder
void fmaxmin(const float *p, long n, float *pmax, float *pmin){
	float max=-INFINITY, min=INFINITY;
	for(long i=0; i<n; i++){
		if(!isnan(p[i])){//not NAN
			if(p[i]>max){
				max=p[i];
			}
			if(p[i]<min){
				min=p[i];
			}
		}
	}
	if(pmax) *pmax=max;
	if(pmin) *pmin=min;
}
/**
   convert float to int
*/
static unsigned int crp(float x, float x0){
	float res=1.5-4.*fabs(x-x0);
	if(res>1) res=1.;
	else if(res<0) res=0.;
	return (unsigned int)(res*255.);
}

/**
   convert float to char with color map*/
void flt2pix(const float *restrict p, void *pix, long nx, long ny, int gray, float *zlim, int zlim_manual, int zlog){
	float max, min;
	fmaxmin(p, nx*ny, &max, &min);
	//info("min=%g, max=%g\n", min, max);
	if(zlog){
		if(max<0){
			float max2=max;
			max=log10(fabs(min));
			min=log10(fabs(max2));
		}else if(min<0){
			max=log10(MAX(fabs(max),fabs(min)));
			min=max-5;
		}else{
			max=log10(max);
			min=log10(min);
		}
		//info("min=%g, max=%g\n", min, max);
	}
	round_limit(&min, &max, 0);
	if(zlim[0]>=zlim[1]||!isfinite(zlim[0])||!isfinite(zlim[1])){//invalid values.
		//info("zlim=[%g,%g] is invalid. replace with %g, %g\n", zlim[0], zlim[1], min, max);
		zlim[0]=min;
		zlim[1]=max;
	}else if(!zlim_manual){
		//update if range change by +15%, -30%
		float thres=(zlim[1]-zlim[0]+max-min)*0.1;
		if(max>zlim[1]+thres||max<zlim[1]-thres*2){
			//info("max updated from %g to %g, thres is %g.\n", zlim[1], max, thres);
			zlim[1]=max;
		}
		if(min>zlim[0]+thres*2||min<zlim[0]-thres){
			//info("min updated from %g to %g, thres is %g.\n", zlim[0], min, thres);
			zlim[0]=min;
		}
	}
	//round_limit(zlim, zlim+1, 0);
	min=zlim[0];
	max=zlim[1];
	if(!gray){/*colored */
		int *pi=(int *)pix;
		float scale, offset;
		if(fabs(max-min)>1.e-4*fabs(min)){
			scale=1./(max-min);
			offset=0;
		} else{
			scale=0;
			offset=0.5;
		}
		for(int i=0; i<nx*ny; i++){
			float xi=p[i];
			if(isnan(xi)){
				pi[i]=0;
			} else{
				if(zlog) xi=log10(fabs(xi));
				float x=(xi-min)*scale+offset;
				pi[i]=255<<24|crp(x, 0.75)<<16|crp(x, 0.5)<<8|crp(x, 0.25);
			}
		}
	} else{/*b/w */
		unsigned char *pc=(unsigned char *)pix;
		float scale=255./(max-min);
		for(int i=0; i<nx*ny; i++){
			float xi=p[i];
			if(isnan(xi)){
				pc[i]=0;
			}else{
				if(zlog) xi=log10(fabs(xi));
				pc[i]=(unsigned char)((xi-min)*scale);
			}
		}
	}
}

void *listen_udp(void *dummy){
	(void)dummy;
	dbg_time("listen_dup listening at socket %d\n", udp_sock);
	char *buf=0;
	size_t bufsize=0;
	int counter=0;
	do{
		counter=udp_recv(&udp_client, &buf, &bufsize);
		info_time("listen_udp: %lu bytes received with counter %d.\n", bufsize, counter);
	} while(counter>0);
	return NULL;
}
/* Initialize */
void drawdata_init(drawdata_t *drawdata){
	drawdata->zoomx=1;
	drawdata->zoomy=1;
	drawdata->square=-1;
	drawdata->legendbox=1;
	drawdata->legendcurve=1;
	drawdata->legendoffx=1;
	drawdata->legendoffy=0;
	drawdata->xylog[0]='n';
	drawdata->xylog[1]='n';
	drawdata->cumulast=-1;/*mark as unknown. */
	pthread_mutex_init(&drawdata->mutex, NULL);
}
/**
 * @brief Free content but keep for reuse.
 * 
 * @param drawdata 
 */
static void drawdata_freecontent(drawdata_t *drawdata){
	/*Only free the input received via fifo from draw.c */
	if(!drawdata->fig) return;
	FREE(drawdata->fig);
	FREE(drawdata->name);
	FREE(drawdata->title);
	FREE(drawdata->xlabel);
	FREE(drawdata->ylabel);
	free_strarr(drawdata->legend, drawdata->nptsmax); drawdata->legend=NULL;
	free_strarr(drawdata->legend_ellipsis, drawdata->nptsmax); drawdata->legend_ellipsis=NULL;
	FREE(drawdata->filename);
	FREE(drawdata->filename_gif);
	FREE(drawdata->p0);
	FREE(drawdata->p1);
	FREE(drawdata->p);

	if(drawdata->npts>0){
		for(int ipts=0; ipts<drawdata->nptsmax; ipts++){
			FREE(drawdata->pts[ipts]);
		}
		FREE(drawdata->pts);
		FREE(drawdata->ptsdim);
	}
	if(drawdata->nstylemax>0){
		FREE(drawdata->style);
	}
	if(drawdata->ncirmax>0){
		FREE(drawdata->cir);
	}
	FREE(drawdata->limit_data);
	FREE(drawdata->limit_cumu);
	//info_time("Free content: %p\n", drawdata);
}

static drawdata_t *HEAD=NULL;
static drawdata_t *drawdata_get(char **fig, char **name, int reset){
	drawdata_t *drawdata=0;
	for(drawdata_t **pp=&HEAD; *pp; ){
		drawdata_t *p=*pp;
		if(__atomic_load_n(&p->recycle,__ATOMIC_SEQ_CST)==2){
			drawdata_freecontent(p);
			if(p->io_time+10<myclockd()){//recycled by GUI. give 10 seconds before free to avoid race condition
				*pp=p->next;
				FREE(p);
			}else{
				pp=&p->next;
			}
		} else {
			if(!__atomic_load_n(&p->recycle,__ATOMIC_SEQ_CST) && p->fig && p->name && !strcmp(p->fig, *fig)&&!strcmp(p->name, *name)){
				drawdata=p;
			}
			pp=&p->next;
		}
	}
	if(!drawdata){
		//create a new node
		drawdata=mycalloc(1, drawdata_t);
		drawdata->next=HEAD;
		HEAD=drawdata;
		drawdata_init(drawdata);
		drawdata->fig=*fig; *fig=0;
		drawdata->name=*name; *name=0;
	} else{
		//reset image, npoints, to default. do not reset memory
		/*while(!drawdata->drawn && drawdata->ready){
			warning_time("Wait for previous data to draw before receiving new data\n");
			mysleep(1);
		}*/
		FREE(*fig);
		FREE(*name);
	}
	if(reset){
		drawdata->update_limit=1;
	}
	if(drawdata->session<session){
		drawdata->session=session;
		drawdata->limit_manual=0;
		drawdata->zlim_manual=0;
		//drawdata->cumu=0;
		drawdata->update_zoom=2;
		drawdata->update_limit=1;
	}
	return drawdata;
}
/**
 * Replace common span of phases by ...
 */
static char** char_ellipsis(char *legends[], int npts){
	if(!legends || npts<1 || !legends[0]){
		return NULL;
	}
	char **legends2=mycalloc(npts, char *);
	for(int ipts=0; ipts<npts; ipts++){
		if(legends[ipts]){
			legends2[ipts]=strdup(legends[ipts]);
		}
	}
	int slen=strlen(legends[0]);
	if(npts<1 || slen<20){
		return legends2;
	}
	
	int cstart=-1;//start of common string
	int clen=0;//length if common string
	for(int is=0; is<slen; is++){
		char c0=legends2[0][is];
		int eq=1;
		//check whether this position is common.
		for(int ip=1; ip<npts; ip++){
			if(c0!=legends2[ip][is]){
				eq=0;break;
			}
		}
		
		if(eq){//common
			if(cstart!=-1){//continuation
				clen++;
			}else{
				cstart=is;//mark the start
				clen=0;
			}
		}
		if(cstart!=-1 && (!eq || is+1==slen)){//end of common or end of str.
			if(clen>4){//turn to ellipsis 
				for(int ip=0; ip<npts; ip++){
					for(int j=cstart; j<cstart+3; j++){
						legends2[ip][j]='.';
					}
					memmove(legends2[ip]+cstart+3, legends2[ip]+cstart+clen+1, strlen(legends2[ip])-cstart-clen);
				}
			}
			is=cstart+2;
			slen=strlen(legends2[0]);
			cstart=-1;
		}
	}
	return legends2;
}
/**
 * Delete pages that are not updated between DRAW_INIT and DRAW_FINAL
 * */
static void drawdata_clear_older(float timclear){
	if(!HEAD) return;
	for(drawdata_t *p=HEAD; p; p=p->next){
		if(p->drawarea && p->io_time<timclear){
			info_time("Request deleting page %s %s\n", p->fig, p->name);
			g_idle_add((GSourceFunc)delete_page, p);
		}
	}
}
//CATCH_TO handles time out.
#define READ_BUF_INT(p) p=*(int*)(buf+boff); boff+=sizeof(int)
#define READ_BUF_N(p, len) memcpy(p, buf+boff,len); boff+=len
#define READ_BUF_STR(p) 
#define CATCH_TO(A,p) if(A) {if(errno==EAGAIN ||errno==EWOULDBLOCK){continue;}else{ close(sock); sock=-1; dbg_time("read " #p " failed %s.", strerror(errno)); break;}}
#define CATCH(A,p) if(A) {close(sock); sock=-1; dbg_time("read " #p " failed %s.", strerror(errno)); break;}
#define STREADINT_TO(p) if(boff<blen){READ_BUF_INT(p);} else CATCH_TO(streadint(sock, &p), p)
#define STREADINT(p) if(boff<blen){READ_BUF_INT(p);} else CATCH(streadint(sock, &p),p)
#define STREAD(p,len) if(boff<blen){READ_BUF_N(p, len);} else CATCH(stread(sock,p,len),p)
#define STREADSTR(p)  if(boff<blen){int len; READ_BUF_INT(len); p=(char*)realloc(p, len); READ_BUF_N(p, len);}else CATCH(streadstr(sock, &p),p)
//read and convert incoming data to float
#define STREADFLT(p,len) if(boff<blen){READ_BUF_N(p, len*byte_float);}else CATCH(stread(sock, p, len*byte_float),p) \
if(byte_float==8){								\
	for(int i=0; i<len; i++){					\
	    ((float*)p)[i]=(float)(((double*)p)[i]);\
	}											\
}
int sock;//socket
int client_pid=-1;//client PID. -1: disconnected. 0: idle. >1: active plotting
int client_pidold=-1; //old client PID.
int draw_id=0; //1: maos, 2: drawres.
int keep_listen=1;//set to 0 to stop listening
int draw_single=0;//whether client only wants to draw to the active tab.
static drawdata_t *drawdata=NULL;//current
static drawdata_t *drawdata_prev=NULL;//previous
char *client_hostname=NULL;
char *client_path=NULL;//received
char *client_path_full=NULL;//with ~ expanded.
char *client_exename=NULL;//client executable name.
int npts=0;
void *listen_draw(void *user_data){
	if(user_data){
		char *str2=0;
		sock=strtol((char *)user_data, &str2, 10);
		if(str2!=user_data){//argument is a number
			if(sock<0){
				error("sock=%d is invalid\n", sock);
			}
			client_hostname=strdup(addr2name(socket_peer(sock)));
		} else{//not a number, hostname
			client_hostname=strdup((char *)user_data);
			sock=-1;
		}
	}//else: use previous information
	while(keep_listen){
		client_pid=-1;
		if(sock<0&&client_hostname){
			dbg_time("Connecting to %s\n", client_hostname);
			sock=scheduler_connect(client_hostname);
			if(sock==-1){
				warning_time("connect to %s failed, retry in 60 seconds.\n", client_hostname);
				mysleep(60);
				continue;
			}
			int cmd[2]={CMD_DRAWSER, -draw_id};
			if(stwriteintarr(sock, cmd, 2)||streadintarr(sock, cmd, 1)||cmd[0]){
				warning_time("Failed to register sock in scheduler, retry in 60 seconds.\n");
				close(sock);
				mysleep(60);
				continue;
			}
		}
		g_idle_add(update_fpslabel, GINT_TO_POINTER(0));//run once
		if(sock!=-1){
			client_pid=0;
			//we set socket timeout to check disconnection.
			//server sends heartbeat every 10 seconds (since 2021-09-29).
			if(socket_block(sock, 0)||socket_recv_timeout(sock, 30)){//was 600. changed to 30 to detect disconnection.
				sock=-1;
				warning("Set sock block and timout failed.\n");
			}
		}
		draw_single=0;
		char *fig=0;
		char *name=0;
		int cmd=0;
		int nlen=0;
		char *buf=0;//read socket into buffer
		int boff=0;//buffer read offset
		int blen=0;//content length
		int bsize=0;//memory allocation size
		//int pid=getpid();
		if(sock!=-1) dbg_time("listen_draw is listening at %d\n", sock);
		while(sock!=-1){
			if(boff>blen){
				error("buffer over use: boff=%d, blen=%d, bsize=%d\n", boff, blen, bsize);
			}
			STREADINT_TO(cmd);//handles time out.
			if(cmd==DRAW_ENTRY){//every message in new format start with DRAW_ENTRY.
				STREADINT(nlen);
				STREADINT(cmd);
				if(1 && nlen && boff>=blen){//read into buffer
					if(nlen>bsize){
						buf=(char*)realloc(buf, nlen);
						bsize=nlen;
					}
					STREAD(buf, nlen);
					blen=nlen;
					boff=0;//reset position in buffer
				}
			}
			//if(cmd!=DRAW_HEARTBEAT) dbg_time("%d received %d\n", pid, cmd);

			switch(cmd){
			case DRAW_FRAME:{//in new format, every frame start with this. Place holder to handle UDP.
				int sizes[4];
				STREAD(sizes, sizeof(int)*4);
			};break;
			case DRAW_START://0
				//tic;
				if(drawdata) warning_time("listen_draw: drawdata=%p, should be NULL.\n", drawdata);
				if(fig) warning_time("fig=%s, should be NULL.\n", fig);
				if(name) warning_time("name=%s\n, should be NULL.\n", name);
				break;
			case DRAW_DATA:/*image data. */
			{
				int32_t header[2];
				STREAD(header, 2*sizeof(int32_t));
				long tot=header[0]*header[1];
				if(drawdata->nmax<tot){
					drawdata->nmax=tot;
					drawdata->p0=(float*)realloc(drawdata->p0, tot*byte_float);//use byte_float to avoid overflow
					if(lpf<1){
						free(drawdata->p1);//recreate below; 
						drawdata->p1=NULL;
					}
					
					drawdata->p=(char*)realloc(drawdata->p, tot*sizeof(int));
					
				}
				if(tot>0){
					if(lpf<1&&drawdata->p1){
						STREADFLT(drawdata->p1, tot);
						for(long i=0; i<tot; i++){
							drawdata->p0[i]=drawdata->p0[i]*(1.-lpf)+drawdata->p1[i]*lpf;
						}
					}else{
						STREADFLT(drawdata->p0, tot);
						if(lpf<1&&!drawdata->p1){
							drawdata->p1=(float*)realloc(drawdata->p1, tot*byte_float);
						}
					}
				}
				drawdata->nx=header[0];
				drawdata->ny=header[1];
			}
			break;
			case DRAW_HEARTBEAT:/*no action*/
				io_heartbeat=myclocki();
				//dbg_time("heatbeat=%d\n", io_heartbeat);
				break;
			case DRAW_POINTS:
			{
				int ipts=npts;
				npts++;
				if(npts>drawdata->nptsmax){
					drawdata->pts=myrealloc(drawdata->pts, npts, float *);
					drawdata->style_pts=myrealloc(drawdata->style_pts, npts, int);
					drawdata->ptsdim=(int(*)[2])realloc(drawdata->ptsdim, npts*2*sizeof(int));
					drawdata->legend=myrealloc(drawdata->legend, npts, char *);
					for(; drawdata->nptsmax<npts; drawdata->nptsmax++){
						drawdata->pts[drawdata->nptsmax]=NULL;
						drawdata->style_pts[drawdata->nptsmax]=0;
						drawdata->ptsdim[drawdata->nptsmax][0]=0;
						drawdata->ptsdim[drawdata->nptsmax][1]=0;
						drawdata->legend[drawdata->nptsmax]=NULL;
					}
				}
				int nptsx, nptsy;

				STREADINT(nptsx);
				STREADINT(nptsy);
				STREADINT(drawdata->square);
				drawdata->grid=1;
				if(drawdata->ptsdim[ipts][0]*drawdata->ptsdim[ipts][1]<nptsx*nptsy){
					drawdata->pts[ipts]=(float*)realloc(drawdata->pts[ipts], nptsx*nptsy*byte_float);
					//drawdata->icumu=0;
				}
				drawdata->ptsdim[ipts][0]=nptsx;
				drawdata->ptsdim[ipts][1]=nptsy;
				if(nptsx*nptsy>0){
					STREADFLT(drawdata->pts[ipts], nptsx*nptsy);

#define FIND_ICUMU(T)\
	T sumc=0; T *p=(T*)drawdata->pts[ipts];\
		for(isim0=nptsx-1; isim0>1; isim0--){\
			sumc+=p[isim0];\
			if((nptsx-isim0)>25 && p[isim0-1]*(nptsx-isim0)>sumc*1.2){\
				break;\
			}\
		}
					if(!drawdata->icumu||drawdata->icumu>nptsx){
						int isim0=0;
						if(byte_float==4){
							FIND_ICUMU(float);
						}else if(byte_float==8){
							FIND_ICUMU(double);
						}else{
							isim0=nptsx/5;
						}
						drawdata->icumu=isim0;
					}
				}
				//info("%s %s: %dx%d\n", drawdata->fig, drawdata->name, nptsx, nptsy);
			}
			break;
			case DRAW_STYLE:
				STREADINT(drawdata->nstyle);
				if(drawdata->nstylemax<drawdata->nstyle){
					drawdata->style=myrealloc(drawdata->style, drawdata->nstyle, uint32_t);
					drawdata->nstylemax=drawdata->nstyle;
				}
				STREAD(drawdata->style, sizeof(uint32_t)*drawdata->nstyle);
				break;
			case DRAW_CIRCLE:
				STREADINT(drawdata->ncir);
				if(drawdata->ncirmax<drawdata->ncir){
					drawdata->cir=(float(*)[4])realloc(drawdata->cir, 4*drawdata->ncir*byte_float);
					drawdata->ncirmax=drawdata->ncir;
				}
				STREADFLT(drawdata->cir, 4*drawdata->ncir);
				break;
			case DRAW_LIMIT:
				if(!drawdata->limit_data){
					drawdata->limit_data=(float*)malloc(4*byte_float);
				}
				STREADFLT(drawdata->limit_data, 4);
				drawdata->limit_manual=1;
				break;
			case DRAW_FIG:
				STREADSTR(fig);
				break;
			case DRAW_NAME:
				STREADSTR(name);
				if(fig&&name){
					drawdata=drawdata_get(&fig, &name, 0);
					npts=0;
				} else{
					warning_time("Invalid usage: fig should be provided before name.\n");
				}
				break;
			case DRAW_TITLE:
				STREADSTR(drawdata->title);
				break;
			case DRAW_XLABEL:
				STREADSTR(drawdata->xlabel);
				break;
			case DRAW_YLABEL:
				STREADSTR(drawdata->ylabel);
				break;
			case DRAW_ZLIM:
				if(drawdata->zlim_manual!=2){
					STREADFLT(drawdata->zlim, 2);
					drawdata->zlim_manual=1;
					drawdata->zlog_last=0;//zlim is not in log format.
				}else{
					double temp[2];
					STREADFLT(temp, 2);
				}
				break;
			case DRAW_LEGEND:
				for(int i=0; i<npts; i++){
					STREADSTR(drawdata->legend[i]);
				}
				break;
			case DRAW_XYLOG:
				STREAD(drawdata->xylog, sizeof(char)*2);
				break;
			case DRAW_FINAL:
				session++;
				//dbg_time("client is done\n");
				if(io_timeclear){
					drawdata_clear_older(io_timeclear);
				}
				client_pid=0;
				g_idle_add(finalize_gif, NULL);
				break;
			case DRAW_FLOAT:
				//notice that this value can change from plot to plot
				//currently, points uses real (default to double), while image uses float.
				//memory in drawdaemon are allocated ALWAYS using float
				STREADINT(byte_float);
				if(byte_float>8){
					error("invalid byte_float=%d\n", byte_float);
				}
				break;
			case DRAW_UDPPORT://received UDP port from client
			{
				STREADINT(client_port);
				client_addr=socket_peer(sock);
				info_time("received udp port %d fron client %s\n", client_port, addr2name(client_addr));

				if(udp_sock<=0){
					udp_sock=bind_socket(SOCK_DGRAM, 0, 0);
				}
				int server_port=socket_port(udp_sock);
				struct sockaddr_in add;
				add.sin_family=AF_INET;
				add.sin_addr.s_addr=client_addr;
				add.sin_port=htons(client_port);
				if(connect(udp_sock, (const struct sockaddr *)&add, sizeof(add))){
					warning_time("connect udp socket to client failed with error %d\n", errno);
				} else{
					//initial handshake with fixed buffer size of 64 ints. The length can not be increased.
					int cmd2[64]={0};
					cmd2[0]=DRAW_ENTRY;
					cmd2[1]=sizeof(int)*4;
					cmd2[2]=1;//protocol version
					cmd2[3]=server_port;
					cmd2[4]=UDP_PAYLOAD;
					cmd2[5]=UDP_HEADER;
					udp_client.header=UDP_HEADER;
					udp_client.payload=UDP_PAYLOAD;
					udp_client.peer_addr=add;
					udp_client.version=1;
					udp_client.sock=udp_sock;
					if(send(udp_sock, cmd2, sizeof(cmd2), 0)<(ssize_t)sizeof(cmd2)){
						warning_time("write to client failed with error %d\n", errno);
					} else{
						thread_new(listen_udp, NULL);
					}
				}
			}
			break;
			case DRAW_INIT:
			{
				io_timeclear=myclockd();//when DRAW_FINAL is called, plots older than this will will be cleared
			}
			break;
			case DRAW_PID:
			{
				STREADINT(client_pid);
				client_pidold=client_pid;
				g_timeout_add(2000, update_fpslabel, GINT_TO_POINTER(1));//will cancel when idle or disconnected
			}
			break;
			case DRAW_ZLOG://flip zlog
				STREADINT(drawdata->zlog);
				break;
			case DRAW_PATH:
				STREADSTR(client_path);
				if(client_path && client_path[0]=='~'){
					client_path_full=stradd(HOME, client_path+1, NULL);
				}
				g_idle_add(update_title, NULL);
				break;
			case DRAW_EXENAME:
				STREADSTR(client_exename);
				break;
			case DRAW_END:
			{
				drawdata->npts=npts;
				if(drawdata->npts>0){
					drawdata->cumuquad=1;
					if(drawdata->nstyle>1){
						if(drawdata->nstyle!=drawdata->npts){
							warning_time("nstyle must equal to npts\n");
							drawdata->nstyle=0;/*disable it. */
							free(drawdata->style);
						}
					}
					drawdata->update_limit=1; //data range may be changed. recompute.
					if(drawdata->legend){
						drawdata->legend_ellipsis=char_ellipsis(drawdata->legend, npts);
					}
				}
				if(drawdata->nx&&drawdata->ny){/*draw image */
					if(!drawdata->limit_data){
						drawdata->limit_data=mycalloc(4, float);
						drawdata->limit_manual=0;
					}
					if(!drawdata->limit_manual){
						drawdata->limit_data[0]=-0.5;
						drawdata->limit_data[1]=drawdata->nx-0.5;
						drawdata->limit_data[2]=-0.5;
						drawdata->limit_data[3]=drawdata->ny-0.5;
					}
					if(drawdata->nx_last!=drawdata->nx || drawdata->ny_last!=drawdata->ny){
						drawdata->update_zoom=2;
						drawdata->nx_last=drawdata->nx;
						drawdata->ny_last=drawdata->ny;
					}
					if(drawdata->square==-1) drawdata->square=1;//default to square for images.
					drawdata->zlim_changed=1;//ask cairo_draw to reconvert the data
				}
				drawdata->frame_io++;
				drawdata->ready=1;
				
				if(drawdata_prev&&drawdata_prev==drawdata){//same drawdata is updated, enable computing framerate.
					io_time2=io_time1;
				} else{
					io_time2=0;//this disables frame rate printing
				}
				io_time1=myclockd();
				drawdata->io_time=io_time1;
				drawdata_prev=drawdata;//for computing time
				g_idle_add((GSourceFunc)addpage, drawdata);//race condition: it may happen after drawdata was destroyed in the GUI.
				//info("Queuing %p: %s, %s\n", drawdata, drawdata->fig, drawdata->name);
				drawdata=NULL;
			}
			break;
			case -1://read failed or server request close.
				close(sock);
				sock=-1;
				break;
			default:
				warning_time("Unknown cmd: %d with size %d\n", cmd, nlen);
				if(nlen && blen){
					boff=blen;//ignore the message
				}
			}/*switch */
			cmd=-1;
			nlen=0;
		}/*while */
	}
	g_idle_add(update_fpslabel, GINT_TO_POINTER(0));//run once
	dbg_time("Stop listening.\n");
	if(sock!=-1) close(sock);
	sock=-1;
	client_pid=-1;
	return NULL;
}
