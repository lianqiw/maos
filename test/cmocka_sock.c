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
/*
  Test socket helpers (sys/sock.c) and socket IO (sys/sockio.c, readsock and
  writesock in math/cell.c). Only pipes, socketpairs, and the loopback
  interface (127.0.0.1 with ephemeral ports) are used. Every test arms an
  alarm so that a blocking call cannot hang the program.
*/
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <poll.h>
#include <arpa/inet.h>
#include <netinet/tcp.h>
#include "unittest.h"

#define ARM_ALARM alarm(30) /*SIGALRM kills the test if any IO blocks for too long*/
#define DISARM_ALARM alarm(0)

/** Create a pair of connected unix stream sockets. */
static void make_socketpair(int sv[2]){
	if(socketpair(AF_UNIX, SOCK_STREAM, 0, sv)){
		fail_msg("socketpair failed: %s", strerror(errno));
	}
}
/** Create a TCP socket listening on 127.0.0.1 at an ephemeral port. Returns the socket and sets *port. */
static int make_listener(int *port){
	int sock=bind_socket(SOCK_STREAM, "127.0.0.1", 0);
	assert_true(sock>=0);
	assert_int_equal(listen(sock, 4), 0);
	*port=socket_port(sock);
	assert_true(*port>0);
	return sock;
}
/** Return a loopback port that nobody listens on. */
static int closed_port(void){
	int sock=bind_socket(SOCK_STREAM, "127.0.0.1", 0);
	int port=socket_port(sock);
	close(sock);
	return port;
}
/** stwrite/stread and the int helpers transfer data over a pipe, report partial reads and errors. */
static void sock_stio_pipe(void **state){
	(void)state;
	ARM_ALARM;
	int fd[2];
	assert_int_equal(pipe(fd), 0);
	assert_int_equal(issock(fd[0]), 0);
	const char msg[]="hello world";
	assert_int_equal(stwrite(fd[1], msg, sizeof(msg)), 0);
	assert_int_equal(stwriteint(fd[1], -123456), 0);
	int arr[5]={1, -2, 3, -4, 5};
	assert_int_equal(stwriteintarr(fd[1], arr, 5), 0);
	assert_int_equal(wrap_io(fd[1], (void *)msg, 3, 1), 3);//mode 1 is write
	char buf[sizeof(msg)]={0};
	assert_int_equal(stread(fd[0], buf, sizeof(msg)), 0);
	assert_string_equal(buf, msg);
	int val=0;
	assert_int_equal(streadint(fd[0], &val), 0);
	assert_int_equal(val, -123456);
	int arr2[5]={0};
	assert_int_equal(streadintarr(fd[0], arr2, 5), 0);
	assert_memory_equal(arr, arr2, sizeof(arr));
	memset(buf, 0, sizeof(buf));
	assert_int_equal(wrap_io(fd[0], buf, 3, 0), 3);//mode 0 is read
	assert_memory_equal(buf, "hel", 3);
	//partial read: only 2 bytes are available before end of file.
	assert_int_equal(stwrite(fd[1], "ab", 2), 0);
	close(fd[1]);
	assert_int_equal(stread(fd[0], &val, sizeof(int)), -2);
	assert_int_equal(wrap_io(fd[0], buf, 4, 'r'), 0);//end of file
	close(fd[0]);
	//invalid descriptor
	assert_int_equal(stread(fd[0], buf, 4), -1);
	assert_int_equal(stwrite(fd[0], buf, 4), -1);
	assert_int_equal(wrap_io(-1, buf, 4, 'w'), -1);
	errno=EINTR;//used to loop forever with an invalid mode
	assert_int_equal(wrap_io(-1, buf, 4, 'x'), -1);
	//convert_ans
	assert_int_equal(convert_ans(10, 10), 0);
	assert_int_equal(convert_ans(3, 10), -2);
	assert_int_equal(convert_ans(-1, 10), -1);
	assert_int_equal(stwrite(fd[0], buf, 0), 0);//nothing to transfer
	DISARM_ALARM;
}
/** Strings and string arrays survive stwritestr/streadstr over a socket pair. */
static void sock_stio_str(void **state){
	(void)state;
	ARM_ALARM;
	int sv[2];
	make_socketpair(sv);
	assert_int_equal(issock(sv[0]), 1);
	assert_int_equal(issock(-1), 0);
	char *str=NULL;
	assert_int_equal(stwritestr(sv[0], "a string"), 0);
	assert_int_equal(streadstr(sv[1], &str), 0);
	assert_string_equal(str, "a string");
	free(str);
	//NULL and empty strings
	str=(char *)"x";
	assert_int_equal(stwritestr(sv[0], NULL), 0);
	assert_int_equal(streadstr(sv[1], &str), 0);
	assert_null(str);
	assert_int_equal(stwritestr(sv[0], ""), 0);
	assert_int_equal(streadstr(sv[1], &str), 0);
	assert_non_null(str);
	assert_string_equal(str, "");
	free(str);
	//string array with an empty entry
	const char *arr[]={"first", "", "third"};
	assert_int_equal(stwritestrarr(sv[0], arr, 3), 0);
	char **arr2=NULL;
	int narr=0;
	assert_int_equal(streadstrarr(sv[1], &arr2, &narr), 0);
	assert_int_equal(narr, 3);
	assert_string_equal(arr2[0], "first");
	assert_string_equal(arr2[1], "");
	assert_string_equal(arr2[2], "third");
	free_strarr(arr2, narr);
	//an empty array
	assert_int_equal(stwritestrarr(sv[0], arr, 0), 0);
	assert_int_equal(streadstrarr(sv[1], &arr2, &narr), 0);
	assert_int_equal(narr, 0);
	free(arr2);
	//negative count is rejected
	assert_int_equal(stwriteint(sv[0], -3), 0);
	assert_int_equal(streadstrarr(sv[1], &arr2, &narr), -1);
	//truncated string: length 10 but only 3 bytes before the connection is closed.
	assert_int_equal(stwriteint(sv[0], 10), 0);
	assert_int_equal(stwrite(sv[0], "abc", 3), 0);
	shutdown(sv[0], SHUT_WR);
	str=NULL;
	assert_int_equal(streadstr(sv[1], &str), -2);
	assert_null(str);
	//reading at end of file fails for both
	assert_int_not_equal(streadstr(sv[1], &str), 0);
	assert_null(str);
	assert_int_not_equal(streadstrarr(sv[1], &arr2, &narr), 0);
	close(sv[1]);
	//writing to a closed peer fails with EPIPE without raising SIGPIPE.
	assert_int_equal(stwriteint(sv[0], 1), -1);
	assert_int_not_equal(stwritestr(sv[0], "abc"), 0);
	assert_int_not_equal(stwritestrarr(sv[0], arr, 3), 0);
	close(sv[0]);
	DISARM_ALARM;
}
/** A file descriptor sent with stwritefd() is usable after streadfd(). */
static void sock_fd_pass(void **state){
	(void)state;
	ARM_ALARM;
	int sv[2];
	make_socketpair(sv);
	int fd[2];
	assert_int_equal(pipe(fd), 0);
	assert_int_equal(stwritefd(sv[0], fd[1]), 0);
	int fd2=-1;
	assert_int_equal(streadfd(sv[1], &fd2), 0);
	assert_true(fd2>=0);
	assert_int_not_equal(fd2, fd[1]);
	close(fd[1]);//the received descriptor keeps the pipe open
	assert_int_equal(stwriteint(fd2, 77), 0);
	close(fd2);
	int val=0;
	assert_int_equal(streadint(fd[0], &val), 0);
	assert_int_equal(val, 77);
	assert_int_equal(stread(fd[0], &val, sizeof(int)), -2);//all writers are closed
	//sendmsg and recvmsg fail on a pipe
	assert_int_equal(stwritefd(fd[0], 0), -1);
	assert_int_equal(streadfd(fd[0], &fd2), -1);
	close(fd[0]);
	close(sv[0]);
	close(sv[1]);
	DISARM_ALARM;
}
/** stcheck() detects closed connections without consuming data. */
static void sock_stcheck(void **state){
	(void)state;
	ARM_ALARM;
	int sv[2];
	make_socketpair(sv);
	assert_int_equal(stcheck(sv[1]), 0);//open, no data
	assert_int_equal(stwriteint(sv[0], 5), 0);
	assert_int_equal(stcheck(sv[1]), 0);//open, has data
	assert_int_equal(stcheck(sv[1]), 0);//data is not consumed
	int val=0;
	assert_int_equal(streadint(sv[1], &val), 0);
	assert_int_equal(val, 5);
	close(sv[0]);
	assert_int_equal(stcheck(sv[1]), 1);//orderly shutdown
	close(sv[1]);
	int fd[2];
	assert_int_equal(pipe(fd), 0);
	assert_int_equal(stcheck(fd[0]), 1);//not a socket
	close(fd[0]);
	close(fd[1]);
	DISARM_ALARM;
}
/** Socket option helpers change the options they claim to and reject invalid descriptors. */
static void sock_options(void **state){
	(void)state;
	int sv[2];
	make_socketpair(sv);
	int sock=sv[0];
	assert_int_equal(socket_block(sock, 0), 0);
	assert_true(fcntl(sock, F_GETFL)&O_NONBLOCK);
	//non-blocking read returns immediately with an error
	int val;
	assert_int_equal(stread(sock, &val, sizeof(int)), -1);
	assert_true(errno==EAGAIN||errno==EWOULDBLOCK);
	assert_int_equal(socket_block(sock, 1), 0);
	assert_false(fcntl(sock, F_GETFL)&O_NONBLOCK);
	assert_int_equal(socket_block(-1, 1), -1);

	assert_int_equal(socket_nopipe(sock), 0);
	assert_int_equal(socket_nopipe(-1), -1);

	assert_int_equal(socket_recv_timeout(sock, 3), 0);
	struct timeval tv={0, 0};
	socklen_t len=sizeof(tv);
	assert_int_equal(getsockopt(sock, SOL_SOCKET, SO_RCVTIMEO, &tv, &len), 0);
	assert_int_equal(tv.tv_sec, 3);
	assert_int_equal(socket_send_timeout(sock, 4), 0);
	len=sizeof(tv);
	assert_int_equal(getsockopt(sock, SOL_SOCKET, SO_SNDTIMEO, &tv, &len), 0);
	assert_int_equal(tv.tv_sec, 4);
	assert_int_equal(socket_recv_timeout(-1, 1), -1);
	assert_int_equal(socket_send_timeout(-1, 1), -1);
	//socket only options fail on a pipe
	int fd[2];
	assert_int_equal(pipe(fd), 0);
	assert_true(socket_recv_timeout(fd[0], 1)<0);
	assert_true(socket_send_timeout(fd[0], 1)<0);
	assert_int_equal(socket_port(fd[0]), -1);
	assert_int_equal(socket_peer(fd[0]), 0);
	close(fd[0]);
	close(fd[1]);
	close(fd[0]);//already closed
	assert_int_equal(socket_block(fd[0], 1), -1);//fcntl fails

	cloexec(sv[1]);
	assert_true(fcntl(sv[1], F_GETFD)&FD_CLOEXEC);
	close(sv[0]);
	close(sv[1]);
}
/** bind_socket, connect_port, socket_port, and socket_peer work on the loopback interface. */
static void sock_tcp(void **state){
	(void)state;
	ARM_ALARM;
	int port=0;
	int ls=make_listener(&port);
	assert_int_equal(socket_peer(ls), 0);//not connected
	//connect by address and port
	int cs=connect_port("127.0.0.1", port, 0, 1);
	assert_true(cs>=0);
	int ss=accept(ls, NULL, NULL);
	assert_true(ss>=0);
	assert_int_equal(socket_peer(ss), htonl(INADDR_LOOPBACK));
	assert_int_equal(socket_peer(cs), htonl(INADDR_LOOPBACK));
	assert_int_equal(socket_port(ss), port);
	assert_string_equal(addr2name(htonl(INADDR_LOOPBACK)), "127.0.0.1");
	int one=0;
	socklen_t len=sizeof(one);
	assert_int_equal(getsockopt(cs, IPPROTO_TCP, TCP_NODELAY, &one, &len), 0);
	assert_int_not_equal(one, 0);
	//the returned socket is blocking
	assert_false(fcntl(cs, F_GETFL)&O_NONBLOCK);
	double data[3]={1.5, -2, 1e300};
	assert_int_equal(stwrite(cs, data, sizeof(data)), 0);
	double data2[3];
	assert_int_equal(stread(ss, data2, sizeof(data2)), 0);
	assert_memory_equal(data, data2, sizeof(data));
	assert_int_equal(stwritestr(ss, "reply"), 0);
	char *str=NULL;
	assert_int_equal(streadstr(cs, &str), 0);
	assert_string_equal(str, "reply");
	free(str);
	close(ss);
	{//peer closed. The FIN may take a moment to arrive (e.g., in macOS).
		struct pollfd pfd={cs, POLLIN, 0};
		poll(&pfd, 1, 5000);
	}
	assert_int_equal(stcheck(cs), 1);
	close(cs);

	//connect with host:port. The port argument is ignored.
	char hostport[64];
	snprintf(hostport, sizeof(hostport), "127.0.0.1:%d", port);
	cs=connect_port(hostport, 1, 0, 0);
	assert_true(cs>=0);
	ss=accept(ls, NULL, NULL);
	assert_true(ss>=0);
	len=sizeof(one);
	assert_int_equal(getsockopt(cs, IPPROTO_TCP, TCP_NODELAY, &one, &len), 0);
	assert_int_equal(one, 0);
	assert_int_equal(stwriteint(ss, 99), 0);
	int val=0;
	assert_int_equal(streadint(cs, &val), 0);
	assert_int_equal(val, 99);
	close(ss);
	close(cs);
	close(ls);

	//non blocking connect to a closed port fails immediately
	int port2=closed_port();
	assert_int_equal(connect_port("127.0.0.1", port2, 0, 0), -1);
	//local unix socket that does not exist
	assert_int_equal(connect_port("/nonexistent-maos-dir/maos.sock", 0, 0, 0), -1);

	//bind to an explicit port enables reuse of address.
	int ls2=bind_socket(SOCK_STREAM, "127.0.0.1", (uint16_t)port2);
	assert_true(ls2>=0);
	assert_int_equal(socket_port(ls2), port2);
	int reuse=0;
	len=sizeof(reuse);
	assert_int_equal(getsockopt(ls2, SOL_SOCKET, SO_REUSEADDR, &reuse, &len), 0);
	assert_int_not_equal(reuse, 0);
	assert_true(fcntl(ls2, F_GETFD)&FD_CLOEXEC);
	close(ls2);
	//UDP on any address
	int us=bind_socket(SOCK_DGRAM, NULL, 0);
	assert_true(us>=0);
	assert_true(socket_port(us)>0);
	close(us);
	DISARM_ALARM;
}
/** Host aliases registered by parse_host() are used by lookup_hostaddr, lookup_hostname, and connect_port. */
static void sock_hosts(void **state){
	(void)state;
	ARM_ALARM;
	int port=0;
	int ls=make_listener(&port);
	int nhost0=nhost;
	parse_host("#comment");
	parse_host("");
	assert_int_equal(nhost, nhost0);
	char line[64];
	//127.0.0.3 is on the loopback interface but nothing listens there; the alias redirects to the listener.
	snprintf(line, sizeof(line), "127.0.0.3=127.0.0.1:%d", port);
	parse_host(line);
	parse_host("maostesthost.example.org=127.0.0.5");
	parse_host("maostesthost2");
	assert_int_equal(nhost, nhost0+3);
	assert_string_equal(hosts[nhost0], "127.0.0.3");
	assert_string_equal(hostshort[nhost0], "127.0.0.3");//numerical address is not shortened
	assert_string_equal(hosts[nhost0+1], "maostesthost.example.org");
	assert_string_equal(hostshort[nhost0+1], "maostesthost");
	assert_string_equal(hostshort[nhost0+2], "maostesthost2");
	assert_string_equal(lookup_hostaddr("127.0.0.3"), line+10);
	assert_string_equal(lookup_hostaddr("maostesthost.example.org"), "127.0.0.5");
	assert_string_equal(lookup_hostaddr("maostesthost2"), "maostesthost2");
	assert_string_equal(lookup_hostaddr("unknownhost"), "unknownhost");
	assert_null(lookup_hostaddr(NULL));
	assert_string_equal(lookup_hostname("127.0.0.5"), "maostesthost.example.org");
	assert_string_equal(lookup_hostname("127.0.0.9"), "127.0.0.9");
	assert_null(lookup_hostname(NULL));
	assert_string_equal(addr2name(inet_addr("127.0.0.5")), "maostesthost.example.org");
#if defined(__linux__) //only linux routes all of 127.0.0.0/8 to the loopback interface by default
	//the direct connection to 127.0.0.3 is refused, then the alias is tried.
	int cs=connect_port("127.0.0.3", closed_port(), 0, 0);
	assert_true(cs>=0);
	int ss=accept(ls, NULL, NULL);
	assert_true(ss>=0);
	assert_int_equal(stwriteint(cs, 1234), 0);
	int val=0;
	assert_int_equal(streadint(ss, &val), 0);
	assert_int_equal(val, 1234);
	close(ss);
	close(cs);
#endif
	close(ls);
	DISARM_ALARM;
}

typedef struct{
	int sock;
	const cell *A;
	int nwrite;
}sockwriter_t;
/** Thread that writes A nwrite times to sock with writesock() and shuts down writing if nwrite>1. */
static void *sock_writer(void *data){
	sockwriter_t *w=(sockwriter_t *)data;
	for(int i=0; i<w->nwrite; i++){
		writesock(w->A, w->sock);
	}
	if(w->nwrite>1) shutdown(w->sock, SHUT_WR);
	return NULL;
}
/** Write A nwrite times from a thread and read it back with readsock(). */
static cell *sock_roundtrip(const_anyarray A, int nwrite){
	int sv[2];
	make_socketpair(sv);
	sockwriter_t w={sv[0], A.c, nwrite};
	pthread_t thread;
	pthread_create(&thread, NULL, sock_writer, &w);
	cell *B=readsock(sv[1]);
	pthread_join(thread, NULL);
	//writesock and readsock duplicate the descriptor and leave the socket open.
	assert_int_equal(stwriteint(sv[1], 3), 0);
	int val=0;
	assert_int_equal(streadint(sv[0], &val), 0);
	assert_int_equal(val, 3);
	close(sv[0]);
	close(sv[1]);
	return B;
}
/** Matrices and cells larger than the socket buffer survive writesock/readsock. */
static void sock_readwrite(void **state){
	(void)state;
	ARM_ALARM;
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *A=dnew(200, 300);
	drandn(A, 1, &rstat);
	A->keywords=strdup("dx=0.5\nname=test");
	//A non-cell array is read until end of file and returned as a cell.
	cell *B=sock_roundtrip(A, 2);
	assert_non_null(B);
	assert_int_equal(NX(B), 2);
	assert_int_equal(NY(B), 1);
	for(int i=0; i<2; i++){
		assert_dmat_near(A, dmat_cast(P(B, i)), 0);
		assert_string_equal(P(B, i)->keywords, A->keywords);
	}
	cellfree(B);

	cell *C=cellnew(2, 2);
	P(C, 0)=(cell *)dref(A);
	P(C, 1)=(cell *)cnew(3, 4);
	crandn((cmat *)P(C, 1), 1, &rstat);
	P(C, 3)=(cell *)lnew(5, 1);
	for(int i=0; i<5; i++) P((lmat *)P(C, 3), i)=i*i-7;
	cell *D=sock_roundtrip(C, 1);//cell of known size is read without waiting for end of file
	assert_non_null(D);
	assert_int_equal(NX(D), 2);
	assert_int_equal(NY(D), 2);
	assert_dmat_near(dmat_cast(P(D, 0)), A, 0);
	assert_float_equal(creal(cdiff((cmat *)P(C, 1), cmat_cast(P(D, 1)))), 0, 0);
	assert_true(!P(D, 2)||PN(P(D, 2))==0);
	lmat *L=lmat_cast(P(D, 3));
	assert_non_null(L);
	for(int i=0; i<5; i++) assert_int_equal(P(L, i), i*i-7);
	cellfree(D);
	cellfree(C);
	dfree(A);

	//reading from a closed socket gives NULL
	int sv[2];
	make_socketpair(sv);
	close(sv[0]);
	assert_null(readsock(sv[1]));
	close(sv[1]);
	DISARM_ALARM;
}

/*State shared between listen_port() callbacks and the client thread.*/
typedef struct{
	int done;     //client has finished
	int listener; //tcp listening socket added by the test
	int ticks;    //number of timeout_fun calls
	int nresp;    //number of tcp responses
	int nhttp;    //number of http responses
	int port;     //tcp port
	int fd;       //client end of the socket pair
	int ans[4];   //results observed by the client
}lp_t;
static lp_t lp;
static int lp_responder(struct pollfd *pf, int flag){
	if(flag==-1) return 0;
	int cmd=0;
	if(streadint(pf->fd, &cmd)) return -1;//connection closed
	atomic_add_fetch((unsigned int *)&lp.nresp, 1);
	return stwriteint(pf->fd, cmd+1)?-1:0;
}
static int lp_http_responder(struct pollfd *pf, int flag){
	if(flag==-1) return 0;
	char buf[64]={0};
	if(recv(pf->fd, buf, sizeof(buf)-1, 0)>0&&!mystrcmp(buf, "GET /")){
		atomic_add_fetch((unsigned int *)&lp.nhttp, 1);
		stwrite(pf->fd, "OK", 2);
	}
	return -1;//close the connection
}
static void lp_timeout(void){
	lp.ticks++;
	if(__atomic_load_n(&lp.done, __ATOMIC_SEQ_CST)&&lp.listener!=-1){
		listen_port_del(lp.listener, 1, "test");//closes the listening socket
		lp.listener=-1;
	}
	if(lp.ticks>2000){//about 20 seconds. something is wrong. quit.
		listen_port_del(-2, 1, "test timeout");
	}
}
static void *lp_client(void *data){
	(void)data;
	//a tcp client
	int cs=connect_port("127.0.0.1", lp.port, 1, 1);
	if(cs>=0){
		int val=0;
		if(!stwriteint(cs, 41)&&!streadint(cs, &val)) lp.ans[0]=val;
		//an http client
		int hs=connect_port("127.0.0.1", lp.port, 1, 0);
		if(hs>=0){
			const char *req="GET / HTTP/1.1\r\n\r\n";
			char buf[4]={0};
			if(!stwrite(hs, req, strlen(req))&&!stread(hs, buf, 2)) lp.ans[1]=!strcmp(buf, "OK");
			lp.ans[2]=stread(hs, buf, 1);//closed by the server
			close(hs);
		}
		close(cs);
	}
	//the client through the socket pair
	int val=0;
	if(!stwriteint(lp.fd, 7)&&!streadint(lp.fd, &val)) lp.ans[3]=val;
	close(lp.fd);
	__atomic_store_n(&lp.done, 1, __ATOMIC_SEQ_CST);
	return NULL;
}
/** listen_port() accepts tcp and http clients, serves an existing descriptor, and returns when all are closed. */
static void sock_listen(void **state){
	(void)state;
	ARM_ALARM;
	memset(&lp, 0, sizeof(lp));
	lp.listener=make_listener(&lp.port);
	int sv[2];
	make_socketpair(sv);
	lp.fd=sv[1];
	listen_port_add(lp.listener, POLLIN, NULL, "test_listen");//accept connections internally
	listen_port_add(lp.listener, POLLIN, NULL, "test_listen");//adding again updates the entry
	pthread_t thread;
	pthread_create(&thread, NULL, lp_client, NULL);
	listen_opt_t opt={
		.localpath="relative.sock",//ignored: must be absolute
		.ipaddr=NULL,
		.port=0,
		.fd=sv[0],
		.responder=lp_responder,
		.http_responder=lp_http_responder,
		.timeout_fun=lp_timeout,
		.timeout_sec=0,//raised to 0.001
		.nodelay=0
	};
	listen_port(opt);//returns after all sockets are closed
	register_signal_handler(dummy_signal_handler);//listen_port installs its own
	pthread_join(thread, NULL);
	assert_true(lp.ticks>0);
	assert_true(lp.ticks<=2000);
	assert_int_equal(lp.ans[0], 42);
	assert_int_equal(lp.ans[1], 1);
	assert_int_equal(lp.ans[2], -2);
	assert_int_equal(lp.ans[3], 8);
	assert_int_equal(lp.nresp, 2);
	assert_int_equal(lp.nhttp, 1);
	assert_int_equal(lp.listener, -1);
	assert_int_equal(fcntl(sv[0], F_GETFD), -1);//closed by listen_port
	listen_port_del(12345, 0, "test");//not found
	DISARM_ALARM;
}

static int lq_ticks=0;
static int lq_responder(struct pollfd *pf, int flag){
	(void)pf;
	if(flag==-1) lq_ticks+=1000;//asked to close at quit
	return 0;
}
static void lq_timeout(void){
	lq_ticks++;
	if(lq_ticks==3) raise(SIGUSR2);//caught by the handler of listen_port. Asks it to quit.
}
static int lr_count=0;
static int lr_responder(struct pollfd *pf, int flag){
	(void)pf;
	if(flag==0) lr_count++;
	return -1;//close the socket
}
/** listen_port() listening on a tcp port quits gracefully upon signal and can be started again afterwards. */
static void sock_listen_quit(void **state){
	(void)state;
	ARM_ALARM;
	int sv[2];
	make_socketpair(sv);
	int port=closed_port();
	listen_opt_t opt={
		.localpath=NULL,
		.ipaddr="127.0.0.1",
		.port=port,
		.fd=sv[0],
		.responder=lq_responder,
		.http_responder=NULL,
		.timeout_fun=lq_timeout,
		.timeout_sec=0.01,
		.nodelay=1
	};
	listen_port(opt);
	register_signal_handler(dummy_signal_handler);
	assert_true(lq_ticks>1000);//responder is notified about quitting
	assert_int_equal(fcntl(sv[0], F_GETFD), -1);//closed by listen_port
	//the port is closed
	assert_int_equal(connect_port("127.0.0.1", port, 0, 0), -1);
	close(sv[1]);
	//listen_port can be started again after quitting.
	make_socketpair(sv);
	int val=1;
	assert_int_equal(stwriteint(sv[1], val), 0);
	listen_opt_t opt2={
		.fd=sv[0],
		.responder=lr_responder,
	};
	listen_port(opt2);//returns after the responder closes the only socket
	register_signal_handler(dummy_signal_handler);
	assert_int_equal(lr_count, 1);
	assert_int_equal(fcntl(sv[0], F_GETFD), -1);//closed by listen_port
	close(sv[1]);
	DISARM_ALARM;
}

#define UDP_RESEND 0xEEEE
#define UDP_DONE   0xEEED
/** Create two UDP sockets on 127.0.0.1 connected to each other. */
static void make_udp_pair(int us[2]){
	struct sockaddr_in addr[2];
	for(int i=0; i<2; i++){
		us[i]=bind_socket(SOCK_DGRAM, "127.0.0.1", 0);
		assert_true(us[i]>=0);
		socklen_t len=sizeof(addr[i]);
		assert_int_equal(getsockname(us[i], (struct sockaddr *)&addr[i], &len), 0);
	}
	for(int i=0; i<2; i++){
		assert_int_equal(connect(us[i], (struct sockaddr *)&addr[1-i], sizeof(addr[1-i])), 0);
	}
}
typedef struct{
	udp_t info;
	char *buf;
	size_t len;
	int counter;
	int ans;
}udp_arg_t;
static void *udp_send_thread(void *data){
	udp_arg_t *arg=(udp_arg_t *)data;
	arg->ans=udp_send(&arg->info, arg->buf, arg->len, arg->counter);
	return NULL;
}
static void *udp_recv_thread(void *data){
	udp_arg_t *arg=(udp_arg_t *)data;
	arg->ans=udp_recv(&arg->info, &arg->buf, &arg->len);
	return NULL;
}
/** Receive one packet and check its header. Return the payload size. */
static int udp_get_packet(int sock, char *pkt, int counter, int ipacket, int npacket){
	ssize_t n=recv(sock, pkt, 2048, 0);
	assert_true(n>=12);
	int hdr[3];
	memcpy(hdr, pkt, sizeof(hdr));
	assert_int_equal(hdr[0], counter);
	assert_int_equal(hdr[1], ipacket);
	assert_int_equal(hdr[2], npacket);
	return (int)n-12;
}
/** udp_send() splits data into packets and resends the ones requested. */
static void sock_udp_send(void **state){
	(void)state;
	ARM_ALARM;
	int us[2];
	make_udp_pair(us);
	const size_t len=3500;
	char *data=(char *)malloc(len);
	for(size_t i=0; i<len; i++) data[i]=(char)(i*7+i/256);
	udp_arg_t arg={{{0}, us[0], 1, 1012, 12}, data, len, 9, -100};
	pthread_t thread;
	pthread_create(&thread, NULL, udp_send_thread, &arg);
	socket_recv_timeout(us[1], 10);
	char pkt[2048];
	char *data2=(char *)calloc(len, 1);
	const int npacket=4;
	for(int ip=0; ip<npacket; ip++){
		int n=udp_get_packet(us[1], pkt, 9, ip, npacket);
		assert_int_equal(n, ip+1<npacket?1000:500);
		memcpy(data2+ip*1000, pkt+12, n);
	}
	assert_memory_equal(data, data2, len);
	int garbage[5]={1, 2, 3, 4, 5};
	assert_int_equal(send(us[1], garbage, sizeof(garbage), 0), sizeof(garbage));
	int resend[5]={UDP_RESEND, 9, npacket, 1, 3};
	assert_int_equal(send(us[1], resend, sizeof(resend), 0), sizeof(resend));
	for(int ip=1; ip<3; ip++){
		int n=udp_get_packet(us[1], pkt, 9, ip, npacket);
		assert_int_equal(n, 1000);
		assert_memory_equal(pkt+12, data+ip*1000, n);
	}
	int done[5]={UDP_DONE, 9, npacket, npacket, 0};
	assert_int_equal(send(us[1], done, sizeof(done), 0), sizeof(done));
	pthread_join(thread, NULL);
	assert_int_equal(arg.ans, 0);
	//invalid input
	assert_int_equal(udp_send(NULL, data, len, 0), -1);
	udp_t bad={{0}, 0, 1, 1012, 12};
	assert_int_equal(udp_send(&bad, data, len, 0), -1);
	free(data2);
	free(data);
	close(us[0]);
	close(us[1]);
	DISARM_ALARM;
}
/** Send packet ip of data with the given header. */
static void udp_put_packet(int sock, const char *data, size_t len, int counter, int ip, int npacket){
	char pkt[1012];
	int hdr[3]={counter, ip, npacket};
	memcpy(pkt, hdr, sizeof(hdr));
	size_t n=(size_t)(ip+1)*1000>len?len-(size_t)ip*1000:1000;
	memcpy(pkt+12, data+ip*1000, n);
	assert_int_equal(send(sock, pkt, n+12, 0), (ssize_t)(n+12));
}
/** udp_recv() reassembles out of order packets, ignores foreign and duplicate ones, and requests missing ones. */
static void sock_udp_recv(void **state){
	(void)state;
	ARM_ALARM;
	int us[2];
	make_udp_pair(us);
	const size_t len=3500;
	const int npacket=4;
	char *data=(char *)malloc(len);
	for(size_t i=0; i<len; i++) data[i]=(char)(i*13+5);
	udp_put_packet(us[0], data, len, 5, 0, npacket);
	udp_put_packet(us[0], data, len, 99, 1, npacket);//another frame: dropped
	udp_put_packet(us[0], data, len, 5, 2, npacket);//packet 1 is skipped
	udp_put_packet(us[0], data, len, 5, 0, npacket);//duplicate
	{//packets with out of range index or a short header are dropped.
		int hdr[3]={5, npacket, npacket};
		assert_int_equal(send(us[0], hdr, sizeof(hdr), 0), sizeof(hdr));
		hdr[1]=-1;
		assert_int_equal(send(us[0], hdr, sizeof(hdr), 0), sizeof(hdr));
		hdr[1]=1;
		assert_int_equal(send(us[0], hdr, 8, 0), 8);
	}
	//packet 3 is not sent until requested
	udp_arg_t arg={{{0}, us[1], 1, 1012, 12}, NULL, 0, 0, -100};
	pthread_t thread;
	pthread_create(&thread, NULL, udp_recv_thread, &arg);
	socket_recv_timeout(us[0], 10);
	const int expect[3][5]={
		{UDP_RESEND, 5, npacket, 1, 2},//immediately when packet 2 follows packet 0
		{UDP_RESEND, 5, npacket, 1, 2},//after timeout: the missing packets
		{UDP_RESEND, 5, npacket, 3, 4},
	};
	int cmd[5]={0};
	for(int i=0; i<3; i++){
		assert_int_equal(recv(us[0], cmd, sizeof(cmd), 0), sizeof(cmd));
		assert_memory_equal(cmd, expect[i], sizeof(cmd));
	}
	udp_put_packet(us[0], data, len, 5, 1, npacket);
	udp_put_packet(us[0], data, len, 5, 3, npacket);
	assert_int_equal(recv(us[0], cmd, sizeof(cmd), 0), sizeof(cmd));
	assert_int_equal(cmd[0], UDP_DONE);
	assert_int_equal(cmd[1], 5);
	assert_int_equal(cmd[2], npacket);
	assert_int_equal(cmd[3], npacket);//valid packets
	assert_int_equal(cmd[4], 1);//duplicate packets
	pthread_join(thread, NULL);
	assert_int_equal(arg.ans, 5);
	assert_int_equal(arg.len, npacket*1012);//buffer allocated by udp_recv
	assert_non_null(arg.buf);
	assert_memory_equal(arg.buf, data, len);
	free(arg.buf);
	//invalid input
	size_t len2=10;
	char *buf2=NULL;
	assert_int_equal(udp_recv(NULL, &buf2, &len2), -1);
	assert_int_equal(len2, 0);
	udp_t bad={{0}, us[1], 2, 1012, 12};//unsupported version
	len2=10;
	assert_int_equal(udp_recv(&bad, &buf2, &len2), -1);
	assert_int_equal(len2, 0);
	assert_null(buf2);
	free(data);
	close(us[0]);
	close(us[1]);
	DISARM_ALARM;
}

UNITTEST_MAIN(
	cmocka_unit_test(sock_stio_pipe),
	cmocka_unit_test(sock_stio_str),
	cmocka_unit_test(sock_fd_pass),
	cmocka_unit_test(sock_stcheck),
	cmocka_unit_test(sock_options),
	cmocka_unit_test(sock_tcp),
	cmocka_unit_test(sock_hosts),
	cmocka_unit_test(sock_readwrite),
	cmocka_unit_test(sock_listen),
	cmocka_unit_test(sock_listen_quit),
	cmocka_unit_test(sock_udp_send),
	cmocka_unit_test(sock_udp_recv),
)
