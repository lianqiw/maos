/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
PNEW2(drawdata_mutex);
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
void flt2pix(long nx, long ny, int color, const float *restrict p, void *pout, float *zlim, int zlog){
	float max, min;
	fmaxmin(p, nx*ny, &max, &min);
	if(zlog){
		max=log10(max);
		min=log10(min);
	}
	//update if range change by 10%
	if(max>zlim[1]*1.1||max<zlim[1]*0.9){
		zlim[1]=max;
	}
	if(min>zlim[0]*1.1||min<zlim[0]*0.9){
		zlim[0]=min;
	}
	round_limit(zlim, zlim+1, 0);
	min=zlim[0];
	max=zlim[1];
	if(color){/*colored */
		int *pi=(int *)pout;
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
				if(zlog) xi=log10(xi);
				float x=(xi-min)*scale+offset;
				pi[i]=255<<24|crp(x, 0.75)<<16|crp(x, 0.5)<<8|crp(x, 0.25);
			}
		}
	} else{/*b/w */
		unsigned char *pc=(unsigned char *)pout;
		float scale=255./(max-min);
		for(int i=0; i<nx*ny; i++){
			float xi=p[i];
			if(isnan(xi)){
				pc[i]=0;
			}else{
				if(zlog) xi=log10(xi);
				pc[i]=(unsigned char)((xi-min)*scale);
			}
		}
	}
}
void *listen_udp(void *dummy){
	(void)dummy;
	dbg("listen_dup listening at socket %d\n", udp_sock);
	void *buf=0;
	size_t bufsize=0;
	int counter=0;
	do{
		counter=udp_recv(&udp_client, &buf, &bufsize);
		info("listen_udp: %lu bytes received with counter %d.\n", bufsize, counter);
	} while(counter>0);
	return NULL;
}
void drawdata_free_input(drawdata_t *drawdata){
	/*Only free the input received via fifo from draw.c */

#define FREE(A) free(A); A=NULL;
	FREE(drawdata->p);
	FREE(drawdata->p0);
	FREE(drawdata->limit_data);
	FREE(drawdata->limit_cumu);
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
	FREE(drawdata->fig);
	FREE(drawdata->name);
	FREE(drawdata->title);
	FREE(drawdata->xlabel);
	FREE(drawdata->ylabel);
	if(drawdata->legend){
		for(int i=0; i<drawdata->nptsmax; i++){
			FREE(drawdata->legend[i]);
		}
		FREE(drawdata->legend);
	}

	FREE(drawdata);
}
static drawdata_t *HEAD=NULL;
drawdata_t *drawdata_get(char **fig, char **name, int reset){
	if(!HEAD){
		HEAD=mycalloc(1, drawdata_t);//dummy head for easy handling
	}
	drawdata_t *drawdata=0;
	drawdata_t *ppriv=HEAD;
	for(drawdata_t *p=ppriv->next; p; ppriv=p, p=p->next){
		if(p->recycle){
			ppriv->next=p->next;
			drawdata_free_input(p);
			p=ppriv;
		} else if(!strcmp(p->fig, *fig)&&!strcmp(p->name, *name)){
			drawdata=p;
		}
	}
	if(!drawdata){
		drawdata=mycalloc(1, drawdata_t);
		drawdata->fig=*fig; *fig=0;
		drawdata->name=*name; *name=0;
		drawdata->zoomx=1;
		drawdata->zoomy=1;
		drawdata->square=-1;
		drawdata->format=(cairo_format_t)0;
		drawdata->gray=0;
		drawdata->ticinside=1;
		drawdata->legendbox=1;
		drawdata->legendcurve=1;
		drawdata->legendoffx=1;
		drawdata->legendoffy=0;
		drawdata->xylog[0]='n';
		drawdata->xylog[1]='n';
		drawdata->cumulast=-1;/*mark as unknown. */
		drawdata->limit_manual=0;
		drawdata->next=HEAD->next;
		HEAD->next=drawdata;
	} else if(reset){
		//reset image, npoints, to default. do not reset memory
		/*while(!drawdata->drawn && drawdata->ready){
			warning_time("Wait for previous data to draw before receiving new data\n");
			mysleep(1);
		}*/
		//drawdata->limit_changed=-1;
		//drawdata->drawn=0;
		free(*fig); *fig=0;
		free(*name); *name=0;
		//drawdata->ready=0;
	}
	return drawdata;
}
/**
 * Delete pages that are not updated between DRAW_INIT and DRAW_FINAL
 * */
static void drawdata_clear_older(float timclear){
	if(!HEAD) return;
	for(drawdata_t *p=HEAD->next; p; p=p->next){
		if(p->io_time<timclear){
			p->delete=1;
			info("Request deleting page %s %s\n", p->fig, p->name);
			gdk_threads_add_idle(delete_page, p);
		}
	}
}
#define CATCH(A,p) if(A) {close(sock); sock=-1; dbg("read " #p " failed %s.", strerror(errno)); break;}
#define STREADINT(p) CATCH(streadint(sock, &p),p)
#define STREAD(p,len) CATCH(stread(sock,p,len),p)
#define STREADSTR(p) ({if(p) {free(p);p=NULL;} CATCH(streadstr(sock, &p),p);})
//read and convert incoming data to float
#define STREADFLT(p,len) CATCH(stread(sock, p, len*byte_float),p) \
    if(byte_float!=4){							\
	for(int i=0; i<len; i++){					\
	    ((float*)p)[i]=(float)(((double*)p)[i]);			\
	}								\
    }									
int sock;//socket 
int client_pid=-1;//client PID. -1: disconnected. 0: idle. >1: active plotting
int keep_listen=1;//set to 0 to stop listening
int draw_single=0;//whether client only wants to draw to the active tab.
drawdata_t *drawdata=NULL;//current 
drawdata_t *drawdata_prev=NULL;//previous
char *client_hostname=NULL;
int npts=0;
void *listen_draw(void *user_data){
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
	while(keep_listen){
		client_pid=-1;
		g_idle_add((GSourceFunc)update_title, NULL);
		if(sock<0&&client_hostname){
			dbg_time("Connecting to %s\n", client_hostname);
			sock=scheduler_connect(client_hostname);
			if(sock==-1){
				warning("connect to %s failed, retry in 60 seconds.\n", client_hostname);
				mysleep(60);
				continue;
			}
			int cmd[2]={CMD_DISPLAY, 0};
			if(stwriteintarr(sock, cmd, 2)||streadintarr(sock, cmd, 1)||cmd[0]){
				warning("Failed to register sock in scheduler, retry in 60 seconds.\n");
				close(sock);
				mysleep(60);
				continue;
			}
		}

		if(sock>=0){
			client_pid=0;
			g_idle_add((GSourceFunc)update_title, NULL);
			//we set socket timeout to check disconnection.
			//server sends heartbeat every 10 seconds (since 2021-09-29).
			if(socket_block(sock, 0)||socket_recv_timeout(sock, 600)){
				sock=-1;
			}
		}
		draw_single=0;
		char *fig=0;
		char *name=0;
		int cmd=0;
		int nlen=0;
		//int pid=getpid();
		if(sock!=-1) dbg("listen_draw is listening at %d\n", sock);
		while(sock!=-1){
			STREADINT(cmd);//will block if no data is available.
			if(cmd==DRAW_ENTRY){//every message in new format start with DRAW_ENTRY.
				STREADINT(nlen);
				STREADINT(cmd);
			}
			//if(cmd!=DRAW_HEARTBEAT) dbg_time("%d received %d\n", pid, cmd);

			switch(cmd){
			case DRAW_FRAME:{//in new format, every frame start with this. Place holder to handle UDP.
				int sizes[4];
				STREAD(sizes, sizeof(int)*4);
			};break;
			case DRAW_START:
				//tic;
				if(drawdata) warning("listen_draw: drawdata=%p, should be NULL.\n", drawdata);
				if(fig) warning("fig=%s, should be NULL.\n", fig);
				if(name) warning("name=%s\n, should be NULL.\n", name);
				break;
			case DRAW_DATA:/*image data. */
			{
				int32_t header[2];
				STREAD(header, 2*sizeof(int32_t));
				long tot=header[0]*header[1];
				if(drawdata->nmax<tot){
					drawdata->p0=realloc(drawdata->p0, tot*byte_float);//use byte_float to avoid overflow
					drawdata->p=realloc(drawdata->p, tot*4);
					drawdata->nmax=tot;
				}
				if(tot>0){
					STREADFLT(drawdata->p0, tot);
				}
				drawdata->nx=header[0];
				drawdata->ny=header[1];
				if(drawdata->square==-1) drawdata->square=1;//default to square for images.
			}
			break;
			case DRAW_HEARTBEAT:/*no action*/
				io_heartbeat=myclocki();
				//dbg("heatbeat=%d\n", io_heartbeat);
				break;
			case DRAW_POINTS:
			{
				int ipts=npts;
				npts++;
				if(npts>drawdata->nptsmax){
					drawdata->pts=myrealloc(drawdata->pts, npts, float *);
					drawdata->style_pts=myrealloc(drawdata->style_pts, npts, int);
					drawdata->ptsdim=realloc(drawdata->ptsdim, npts*sizeof(int[2]));
					drawdata->legend=realloc(drawdata->legend, npts*sizeof(char *));
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
					drawdata->pts[ipts]=realloc(drawdata->pts[ipts], nptsx*nptsy*byte_float);
					//drawdata->icumu=0;
				}
				drawdata->ptsdim[ipts][0]=nptsx;
				drawdata->ptsdim[ipts][1]=nptsy;
				if(nptsx*nptsy>0){
					STREADFLT(drawdata->pts[ipts], nptsx*nptsy);
					if(nptsx>50){
						if(!drawdata->icumu){
							drawdata->icumu=nptsx/10;
						}
					}
				}
			}
			break;
			case DRAW_STYLE:
				STREADINT(drawdata->nstyle);
				if(drawdata->nstylemax<drawdata->nstyle){
					drawdata->style=myrealloc(drawdata->style, drawdata->nstyle, int32_t);
					drawdata->nstylemax=drawdata->nstyle;
				}
				STREAD(drawdata->style, sizeof(int32_t)*drawdata->nstyle);
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
					drawdata->limit_data=malloc(4*byte_float);
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
					drawdata=drawdata_get(&fig, &name, 1);
					npts=0;
				} else{
					warning("Invalid usage: fig should be provided before namen");
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
				STREADFLT(drawdata->zlim, 2);
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
				//dbg("client is done\n");
				if(io_timeclear){
					drawdata_clear_older(io_timeclear);
				}
				client_pid=0;
				g_idle_add((GSourceFunc)update_title, NULL);
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
				info("received udp port %d fron client %s\n", client_port, addr2name(client_addr));

				if(udp_sock<=0){
					udp_sock=bind_socket(SOCK_DGRAM, 0, 0);
				}
				int server_port=socket_port(udp_sock);
				struct sockaddr_in add;
				add.sin_family=AF_INET;
				add.sin_addr.s_addr=client_addr;
				add.sin_port=htons(client_port);
				if(connect(udp_sock, (const struct sockaddr *)&add, sizeof(add))){
					warning("connect udp socket to client failed with error %d\n", errno);
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
						warning("write to client failed with error %d\n", errno);
					} else{
						thread_new(listen_udp, NULL);
					}
				}
			}
			break;
			case DRAW_INIT:
			{
				io_timeclear=myclockd();
			}
			break;
			case DRAW_PID:
			{
				STREADINT(client_pid);
				g_idle_add((GSourceFunc)update_title, NULL);
			}
			break;
			case DRAW_END:
			{
				drawdata->npts=npts;
				if(drawdata->npts>0){
					drawdata->cumuquad=1;
					if(drawdata->nstyle>1){
						if(drawdata->nstyle!=drawdata->npts){
							warning("nstyle must equal to npts\n");
							drawdata->nstyle=0;/*disable it. */
							free(drawdata->style);
						}
					}
				}
				if(!drawdata->fig) drawdata->fig=strdup("unknown");
				drawdata->drawn=0;
				drawdata->ready=1;

				if(drawdata_prev&&drawdata_prev==drawdata){//same drawdata is updated, enable computing framerate.
					io_time2=io_time1;
				} else{
					io_time2=0;//this disables frame rate printing
				}
				io_time1=myclockd();
				drawdata->io_time=io_time1;
				drawdata_prev=drawdata;//for computing time
				gdk_threads_add_idle(addpage, drawdata);
				drawdata=NULL;
			}
			break;
			case -1://read failed.
				close(sock); 
				sock=-1;
				break;
			default:
				warning_time("Unknown cmd: %d with size %d\n", cmd, nlen);
				if(nlen){
					void *p=malloc(nlen);
					STREAD(p, nlen);
					free(p);
				}
			}/*switch */
			cmd=-1;
			nlen=0;
		}/*while */
	}
	free(client_hostname);client_hostname=NULL;
	dbg_time("Stop listening.\n");
	if(sock!=-1) close(sock);
	sock=-1;
	client_pid=-1;
	g_idle_add((GSourceFunc)update_title, NULL);
	return NULL;
}
