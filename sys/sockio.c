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
#include <errno.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include <poll.h>
#include "misc.h"
#include "common.h"
#include "sockio.h"
#include "sock.h"
#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#define NEED_NOSIGPIPE 1 /*e.g., macOS: set SO_NOSIGPIPE on the socket instead*/
#endif
/**
   Wrap io function and handle partial operations. Uses send()/recv() for
   sockets and write()/read() otherwise. Retries on EINTR and stops at end of
   file or closed connection.

   @param[in] sfd Socket or file descriptor.
   @param[in,out] p Buffer to write from or read into.
   @param[in] len Number of bytes to transfer.
   @param[in] mode 1 or 'w' for write, 0 or 'r' for read.
   @return Number of bytes successfully transferred, or -1 on error or invalid mode.
*/
ssize_t wrap_io(int sfd, void* p, const size_t len, int mode){
	const int is_write=(mode==1 || mode=='w');
	if(!is_write && !(mode==0 || mode=='r')){
		error("invalid mode=%d (%c)\n", mode, mode);
		return -1;
	}
	int is_sock=issock(sfd);
#if NEED_NOSIGPIPE
	if(is_sock&&is_write){//a closed peer would otherwise raise SIGPIPE and kill the process
		socket_nopipe(sfd);
	}
#endif
	size_t total=0;
	while(total<len){
		ssize_t n=-1;
		if(is_sock){
			if(is_write){
				n=send(sfd, (char*)p+total, len-total, MSG_NOSIGNAL);
			}else{
				n=recv(sfd, (char*)p+total, len-total, 0);//MSG_WAITALL);
			}
		} else {
			if(is_write){
				n=write(sfd, (char*)p+total, len-total);
			}else{
				n=read(sfd, (char*)p+total, len-total);
			}
		}
		if(n==0){
			break;
		}else if(n<0){
			if(errno==EINTR) continue;//retry. do not retry for EGAIN or EWOULDBLOCK
			return -1;
		}else{
			total+=(size_t)n;
		}
	}
	return total;
}
/**
   Write len bytes from p to socket or file sfd.

   @param[in] sfd Socket or file descriptor.
   @param[in] p Data to write.
   @param[in] len Number of bytes to write.
   @return 0 on success, -2 for partial write, and -1 on error.
*/
int stwrite(int sfd, const void* p, const size_t len){
	return convert_ans(wrap_io(sfd, (void*)p, len, 'w'), len);
}
/**
   Read len bytes from socket or file sfd into p.

   @param[in] sfd Socket or file descriptor.
   @param[out] p Buffer to receive the data.
   @param[in] len Number of bytes to read.
   @return 0 on success, -2 for partial read (e.g., connection closed), and -1 on error.
*/
int stread(int sfd, void* p, const size_t len){
	return convert_ans(wrap_io(sfd, p, len, 'r'), len);
}

/**
   Write a string to socket, preceded by its length including the terminating
   NUL. A NULL string is written as length 0.

   @param[in] sfd Socket or file descriptor.
   @param[in] str String to write. May be NULL.
   @return 0 on success, nonzero on failure.
*/
int stwritestr(int sfd, const char* str){
	if(str){
		int len=(int)strlen(str)+1;//notice: if pad to multiple of 4 or 8, must change strlen()+1 in draw.c
		return stwriteint(sfd, len)||stwrite(sfd, str, len);
	} else{
		return stwriteint(sfd, 0);
	}
}
/**
   Read a string written by stwritestr(). The string is always NUL terminated.

   @param[in] sfd Socket or file descriptor.
   @param[out] str Newly allocated string; the caller frees it. Set to NULL for an empty string or on failure.
   @return 0 on success, nonzero on failure, including length above 1 MiB.
*/
int streadstr(int sfd, char** str){
	int len;
	int ans=streadint(sfd, &len);
	if(!ans&&len>(1<<20)){//protect against corrupted or malicious length
		warning("string length %d is too large\n", len);
		ans=-1;
	}
	if(!ans&&len>0){
		*str=mycalloc(len+1, char);//+1 in case the sender did not include NUL
		ans=stread(sfd, *str, (size_t)len);
		if(ans){
			free(*str);
			*str=NULL;
		}
	} else{
		*str=NULL;
	}
	return ans;
}
/**
   Write a string array to socket, preceded by the number of strings.

   @param[in] sfd Socket or file descriptor.
   @param[in] str Array of strings.
   @param[in] nstr Number of strings.
   @return 0 on success, nonzero on failure.
*/
int stwritestrarr(int sfd, const char* const* str, int nstr){
	int ans=stwriteint(sfd, nstr);
	for(int i=0; i<nstr&&!ans; i++){
		ans=stwritestr(sfd, str[i]);
	}
	return ans;
}
/**
   Read a string array written by stwritestrarr().

   @param[in] sfd Socket or file descriptor.
   @param[out] str Newly allocated array of newly allocated strings; the caller frees each string and the array.
   @param[out] nstr Number of strings.
   @return 0 on success, nonzero on failure.
*/
int streadstrarr(int sfd, char*** str, int* nstr){
	int ans=streadint(sfd, nstr);
	if(ans) return ans;
	if(*nstr<0) return -1;
	*str=mycalloc((size_t)*nstr, char*);
	for(int istr=0; istr<*nstr&&!ans; istr++){
		ans=streadstr(sfd, &(*str)[istr]);
	}
	return ans;
}
/**
   Wrie a file descriptor to socket. Transfer a open file (socket) descriptor to
   another process in the same hosts. sfd has to be a AF_UNIX socket.

   @param[in] sfd AF_UNIX socket.
   @param[in] fd File descriptor to transfer.
   @return 0 on success, -1 on failure.
*/
int stwritefd(int sfd, int fd){
	int ans=0;
	char buf[1]={0};
	struct iovec iov;
	iov.iov_base=buf;
	iov.iov_len=1;
	struct msghdr msg={0};
	msg.msg_iov=&iov;
	msg.msg_iovlen=1;
	char cms[CMSG_SPACE(sizeof(int))];
	msg.msg_control=cms;
	msg.msg_controllen=CMSG_LEN(sizeof(int));
	struct cmsghdr* cmsg=CMSG_FIRSTHDR(&msg);
	cmsg->cmsg_level=SOL_SOCKET;
	cmsg->cmsg_type=SCM_RIGHTS;
	cmsg->cmsg_len=CMSG_LEN(sizeof(int));
	memcpy(CMSG_DATA(cmsg), &fd, sizeof(int));
	if(sendmsg(sfd, &msg, 0)<0){
		warning("sendmsg failed to send fd %d over %d\n", fd, sfd);
		ans=-1;
	} else{
		ans=0;
	}
	return ans;
}
/**
   Read a file descriptor from socket. Use with stwritefd();

   @param[in] sfd AF_UNIX socket.
   @param[out] fd Received file descriptor.
   @return 0 on success, -1 on failure.
*/
int streadfd(int sfd, int* fd){
	int ans=0;
	char buf[1]={0};
	struct iovec iov;
	iov.iov_base=buf;
	iov.iov_len=1;
	struct msghdr msg={0};
	msg.msg_iov=&iov;
	msg.msg_iovlen=1;
	char cms[CMSG_SPACE(sizeof(int))];
	msg.msg_control=cms;
	msg.msg_controllen=CMSG_LEN(sizeof(int));
	if(recvmsg(sfd, &msg, 0)<=0){
		warning("recvmsg failed to receive fd over %d\n", sfd);
		ans=-1;
	} else{
		//cmsg is NULL if the peer did not send a descriptor (e.g., over TCP).
		struct cmsghdr* cmsg=CMSG_FIRSTHDR(&msg);
		if(!cmsg || cmsg->cmsg_level!=SOL_SOCKET || cmsg->cmsg_type!=SCM_RIGHTS
			|| cmsg->cmsg_len<CMSG_LEN(sizeof(int))){
			warning("no file descriptor is received over %d\n", sfd);
			ans=-1;
		}else{
			memmove(fd, CMSG_DATA(cmsg), sizeof(int));
			ans=0;
		}
	}
	return ans;
}
/**
   Check whether the client has shutdown a connection (TCP) without blocking
   or consuming data.

   @param[in] sfd Socket.
   @return 1 if closed or failed, 0 otherwise.
*/
int stcheck(int sfd){
	int ans=0;
	char buf[4];
	ssize_t len=recv(sfd, buf, sizeof buf, MSG_DONTWAIT|MSG_PEEK);
	if(len==0){
		dbg_time("Client %d has orderly shutdown\n", sfd);
		ans=1;
	}else if(len==-1 && (errno!=EAGAIN && errno!=EWOULDBLOCK)){
		dbg_time("recv from %d failed with errno %d: %s\n", sfd, errno, strerror(errno));
		ans=1;
	}
	return ans;
}
#define UDP_RESEND 0xEEEE
#define UDP_DONE   0xEEED
/**
   Use UDP socket to send data reliably. Data is split into packets of
   info->payload bytes including header. Afterwards listens for resend
   requests from udp_recv() until it reports completion or times out.

   @param[in] info UDP connection information.
   @param[in] buf Data to send.
   @param[in] len Number of bytes to send.
   @param[in] counter Frame counter written in each packet header.
   @return 0 if sent, -1 if info is invalid.
*/
int udp_send(udp_t* info, char* buf, size_t len, int counter){
	if(!info || info->sock<=0){
		warning("udp_send: info or info->sock is invalid\n");
		return -1;
	}
	TIC;tic;
	dbg("udp_send: %lu bytes with counter %d\n", len, counter);
	size_t header_size=(size_t)info->header;
	size_t packet_size=(size_t)info->payload-header_size;
	size_t npacket=(len+packet_size-1)/packet_size;
	struct msghdr *msg=mycalloc(npacket, struct msghdr);
	struct iovec *iovec=mycalloc(npacket*2, struct iovec);
	char* header=(char*)malloc(npacket*header_size);
	char* buf2=buf;
	char* header2=header;
	ssize_t nsend;
	for(size_t ip=0; ip<npacket; ip++){
		//first iovec contains header for each segment
		iovec[ip*2].iov_base=header2;
		if(info->version==1 && info->header==12){
			((int*)header2)[0]=counter;//frame counter
			((int*)header2)[1]=(int)ip;//subframe counter
			((int*)header2)[2]=(int)npacket;//number of subframes
		}else{
			error("Please implement header information for new version %d\n", info->version);
		}
		header2+=header_size;
		iovec[ip*2].iov_len=header_size;
		iovec[ip*2+1].iov_base=buf2;buf2+=packet_size;
		iovec[ip*2+1].iov_len=packet_size;
		if(ip+1==npacket){
			iovec[ip*2+1].iov_len=len-packet_size*(npacket-1);//last packet may be partial.
		}

		msg[ip].msg_iov=iovec+ip*2;
		msg[ip].msg_iovlen=2;
		size_t msize=iovec[ip*2+1].iov_len+header_size;
		if((nsend=sendmsg(info->sock, msg+ip, 0))<(ssize_t)msize){
			warning("sendmsg failed: %zd bytes sent, expect %zu\n", nsend, msize);
		}
		//usleep(2);
	}
	socket_recv_timeout(info->sock, 10);
	int cmd[5];
	size_t nresent=0, nresent2=0;
	//listen for resend request
	while((nsend=recv(info->sock, cmd, sizeof(cmd), 0))==(ssize_t)sizeof(cmd)){
		if(cmd[0]==UDP_RESEND&&cmd[1]==counter&& cmd[2]==(int)npacket&&cmd[3]>=0&&cmd[4]<=(int)npacket){
			nresent++;
			nresent2+=(size_t)(cmd[4]-cmd[3]);
			for(int ip=cmd[3]; ip<cmd[4]; ip++){
				size_t msize=iovec[ip*2+1].iov_len+header_size;
				if((nsend=sendmsg(info->sock, msg+ip, 0))<(ssize_t)msize){
					warning("sendmsg resend failed: %zd bytes sent, expect %zu\n", nsend, msize);
				}
			}
		} else if(cmd[0]==UDP_DONE&&cmd[1]==counter&&cmd[2]==(int)npacket){
			dbg("udp_send: %d packets of %zu successfully received, %ld total sent, %d received.\n", cmd[3], npacket, npacket+nresent, cmd[3]+cmd[4]);
			break;
		}else{
			dbg("udp_send: garbage received: %d %d %d %d %d\n", cmd[0], cmd[1], cmd[2], cmd[3], cmd[4]);
		}
	}
	dbg("udp_send: %zu resent request received for %zu packets\n", nresent, nresent2);
	//todo: check delivery status.
	free(msg);
	free(iovec);
	free(header);
	toc("udp_send");
	return 0;
}

/**
   Use UDP socket to receive data reliably from udp_send(). Requests resend of
   missing packets and acknowledges when done.

   @param[in] info UDP connection information.
   @param[in,out] pbuf Receive buffer. Reallocated if NULL or too small; the caller frees it.
   @param[in,out] len Size of *pbuf in bytes. Updated if reallocated; set to 0 on failure.
   @return Frame counter of the received data, or -1 on failure.
*/
int udp_recv(udp_t* info, char** pbuf, size_t *len){
	if(!info||info->sock<=0){
		warning("udp_recv info or info->sock is invalid\n");
		*len=0;
		return -1;
	}
	size_t header_size=(size_t)info->header;
	size_t packet_size=(size_t)(info->payload)-header_size;
	//int npacket=(len+packet_size-1)/packet_size;
	size_t npacket;//number of packages
	int counter;//frame counter
	socket_recv_timeout(info->sock, 10);
	ssize_t ans;
	if(info->version==1){
		int meta[3];
		if((ans=recv(info->sock, meta, sizeof(meta), MSG_PEEK))>=(ssize_t)sizeof(meta) && meta[2]>0){//peek at available data.
			info("recv peek got %zd bytes\n", ans);
			counter=meta[0];
			npacket=(size_t)meta[2];
			dbg("%zu packets to be received with counter %d\n", npacket, counter);
		}else{
			warning("Unable to peek data.\n");
			*len=0;
			return -1;
		}
	}else{
		warning("Please implement header information for new version %d\n", info->version);
		*len=0;
		return -1;
	}
	TIC;tic;
	if(!*pbuf||*len<npacket*(size_t)(info->payload)){
		*len=npacket*(size_t)(info->payload);
		warning("Buffer is too small, reallocate to %lu bytes\n", *len);
		if(!(*pbuf=(char*)realloc(*pbuf, *len))){
			warning("realloc buf failed\n");
			*len=0;
			return -1;
		}
	}
	
	struct msghdr* msg=mycalloc(npacket, struct msghdr);
	struct iovec* iovec=mycalloc(npacket*2, struct iovec);
	char* header=(char*)calloc(npacket, header_size);
	char* buf2=*pbuf;
	char* header2=header;
	for(size_t ip=0; ip<npacket; ip++){
		//first iovec contains header for each segment
		iovec[ip*2].iov_base=header2;
		header2+=header_size;
		iovec[ip*2].iov_len=header_size;
		iovec[ip*2+1].iov_base=buf2;buf2+=packet_size;
		iovec[ip*2+1].iov_len=packet_size;
		if(ip+1==npacket){
			iovec[ip*2+1].iov_len=*len-packet_size*(npacket-1);//last packet may be partial.
		}

		msg[ip].msg_iov=iovec+ip*2;
		msg[ip].msg_iovlen=2;
	}
	socket_recv_timeout(info->sock, 1);
	int meta[3];
	size_t nvalid=0;
	int nduplicate=0;//duplicate
	int nresent=0;//resent packet count
	int nresent2=0;//resent request count
	int nreorder=0; //out of order
	int lastind=-1;
	int nretry=0;
retry:	
	nretry++;
	while((ans=recv(info->sock, meta, sizeof(meta), MSG_PEEK))!=-1){//peek at available data.
		int counter2=meta[0];
		int ipacket2=meta[1];
		int npacket2=meta[2];
		if(ans<(ssize_t)sizeof(meta) || counter2!=counter || (size_t)npacket2!=npacket
			|| ipacket2<0 || (size_t)ipacket2>=npacket){
			//drop the packet
			dbg("Drop invalid package %d of %d with counter %d.\n", ipacket2, npacket2, counter2);
			recv(info->sock, meta, sizeof(meta), 0);
		}else{
			int counter3=((int*)msg[ipacket2].msg_iov[0].iov_base)[0];//whether already written
			if(recvmsg(info->sock, msg+ipacket2, 0)!=-1){
				if(counter3){
					nduplicate++;
				}else{//new packet
					nvalid++;
					if(ipacket2<lastind){
						nreorder++;
					}else{
						if(ipacket2>lastind+1){//skipped index
							//dbg("got segment %d after %d. total %d\n", ipacket2, lastind, npacket2);
							int cmd[5]={UDP_RESEND, counter2, npacket2, lastind+1, ipacket2};
							send(info->sock, cmd, sizeof(cmd), 0);
							nresent+=(ipacket2-lastind-1);
						}
						lastind=ipacket2;
					}
				}
			}else{
				dbg("recvmsg failed with %d\n", errno);
			}
		}
	}
	if(nvalid<npacket && nretry<5){
		dbg("%zu packets still missing\n", npacket-nvalid);
		
		int i1=-1;
		for(size_t ip=0; ip<npacket; ip++){
			int counter3=((int*)msg[ip].msg_iov[0].iov_base)[0];
			if(!counter3){
				nresent++;
				if(i1==-1){//mark start of missing date
					i1=(int)ip;
				}
			}else{
				if(i1!=-1){
					int cmd[5]={UDP_RESEND, counter, (int)npacket, i1, (int)ip};//end is exclusive
					send(info->sock, cmd, sizeof(cmd), 0);
					nresent2++;
					i1=-1;
				}
			}
		}
		if(i1!=-1){
			int cmd[5]={UDP_RESEND, counter, (int)npacket, i1, (int)npacket};//end is exclusive
			send(info->sock, cmd, sizeof(cmd), 0);
			nresent2++;
		}
		
		goto retry;
	}
	int cmd[]={UDP_DONE, counter, (int)npacket, (int)nvalid, nduplicate};
	send(info->sock, cmd, sizeof(cmd), 0);
	dbg("%zu of %zu packets received with counter %d. %d duplicated, %d resent, %d reordered.\n", nvalid, npacket, counter, nduplicate, nresent, nreorder);
	dbg("%d resent requestes sent\n", nresent2);
	free(msg);
	free(header);
	free(iovec);
	toc("udp_recv");
	return counter;
}
