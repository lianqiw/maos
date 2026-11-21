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
  Test system utilities in sys/: search path (path.c), process information
  (process.c, process_linux.c), hashing (hashlittle.c, sha1.c), the double
  linked list (dlist.c), zstd streams (zstdio.c), memory management hooks
  (mem.c), and thread helpers (thread.c). Only existing paths are inspected and
  data streams go through pipes.
*/
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <limits.h>
#include "unittest.h"
#include "../sys/zstdio.h"

/** Directories are added to the search path by priority and files are found in them. */
static void sysutil_path(void **state){
	(void)state;
	freepath();
	assert_null(search_file("Makefile.am", 0));
	addpath2(1, "%s", SRCDIR);
	addpath2(2, "%s/sys", SRCDIR);
	addpath("/nonexistent-maos-dir");//ignored
	printpath();
	char *expect_sys=stradd(SRCDIR, "/sys/Makefile.am", NULL);
	char *expect_top=stradd(SRCDIR, "/Makefile.am", NULL);
	//found in both. The higher priority wins.
	char *fn=search_file("Makefile.am", 0);
	assert_non_null(fn);
	assert_string_equal(fn, expect_sys);
	free(fn);
	//moved to the front with higher priority
	addpath2(3, "%s", SRCDIR);
	fn=search_file("Makefile.am", 0);
	assert_string_equal(fn, expect_top);
	free(fn);
	//adding again with lower priority keeps the position
	addpath(SRCDIR);
	fn=search_file("Makefile.am", 0);
	assert_string_equal(fn, expect_top);
	free(fn);
	//relative path below a search directory
	fn=search_file("sys/misc.c", 0);
	char *expect=stradd(SRCDIR, "/sys/misc.c", NULL);
	assert_string_equal(fn, expect);
	free(fn);
	//only in one directory
	fn=search_file("misc.c", 0);
	assert_string_equal(fn, expect);
	free(fn);
	//absolute path
	fn=search_file(expect, 0);
	assert_string_equal(fn, expect);
	free(fn);
	assert_null(search_file("/nonexistent-maos-dir/file", 0));
	assert_null(search_file("nonexistent-maos-file.c", 0));
	assert_null(search_file(NULL, 1));
	//the current directory
	char *cwd=getcwd(NULL, 0);
	fn=search_file(".", 1);
	assert_string_equal(fn, cwd);
	free(fn);
	//remove from path
	rmpath(SRCDIR "/sys");
	rmpath("/nonexistent-maos-dir");//no effect
	assert_null(search_file("misc.c", 0));
	fn=search_file("Makefile.am", 0);
	assert_string_equal(fn, expect_top);
	free(fn);
	//path with trailing / is stored as given
	addpath("/");
	fn=search_file("dev/null", 0);
	assert_string_equal(fn, "/dev/null");
	free(fn);
	freepath();
	assert_null(search_file("Makefile.am", 0));
	free(cwd);
	free(expect);
	free(expect_sys);
	free(expect_top);

	char *config=find_config("maos");
	assert_non_null(config);
	assert_true(isdir(config));
	assert_true(check_suffix(config, "/maos"));
	char *defconf=stradd(config, "/default.conf", NULL);
	assert_true(isfile(defconf));
	free(defconf);
	free(config);
	assert_null(find_config("nonexistent-maos-config"));
	//MAOS_CONFIG_PATH may point to the folder containing config or config itself
	char *env0=getenv("MAOS_CONFIG_PATH");
	if(env0) env0=strdup(env0);
	setenv("MAOS_CONFIG_PATH", SRCDIR "/config", 1);
	config=find_config("maos");
	assert_string_equal(config, SRCDIR "/config/maos");
	free(config);
	assert_null(find_config("nonexistent-maos-config"));
	if(env0){
		setenv("MAOS_CONFIG_PATH", env0, 1);
		free(env0);
	}else{
		unsetenv("MAOS_CONFIG_PATH");
	}
}
#if defined(__linux__)
/** Read the first number from a file in /proc. */
static double read_proc_num(const char *fn){
	double val=-1;
	FILE *fp=fopen(fn, "r");
	if(fp){
		if(fscanf(fp, "%lf", &val)!=1) val=-1;
		fclose(fp);
	}
	return val;
}
#endif
/** Process and system information has sane values. */
static void sysutil_process(void **state){
	(void)state;
	assert_true(NCPU>0);
	assert_int_equal(get_ncpu(), sysconf(_SC_NPROCESSORS_ONLN));
	assert_true(MAXTHREAD>0&&MAXTHREAD<=get_ncpu());
	assert_true(NTHREAD>0);
	assert_int_equal(TCK, sysconf(_SC_CLK_TCK));
	assert_int_equal(PID, getpid());
	assert_true(NMEM>0);
	assert_non_null(HOME);
	assert_non_null(USER);
	assert_non_null(HOST);
	assert_null(strchr(HOST, '.'));//short host name
	assert_non_null(TEMP);
	assert_non_null(DIRSTART);
	assert_non_null(DIREXE);
	assert_non_null(EXENAME);
	assert_string_equal(EXENAME, "cmocka_sysutil");

	char exe[PATH_MAX]={0};
	char link[PATH_MAX]={0};//full path of this executable
	assert_int_equal(get_job_progname(link, PATH_MAX, 0), 0);
#if defined(__linux__)
	ssize_t nlink=readlink("/proc/self/exe", exe, PATH_MAX-1);
	assert_true(nlink>0);
	assert_string_equal(exe, link);
#endif
	strcpy(exe, link);
	assert_true(check_suffix(exe, "/cmocka_sysutil"));
	char *dir=mydirname(exe);
	assert_string_equal(dir, DIREXE);
	free(dir);
	memset(exe, 0, sizeof(exe));
	assert_int_equal(get_job_progname(exe, PATH_MAX, getpid()), 0);
	assert_string_equal(exe, link);
	char shortbuf[8];
	assert_int_equal(get_job_progname(shortbuf, sizeof(shortbuf), 0), 0);//truncated
	assert_int_equal(strlen(shortbuf), 7);
	assert_memory_equal(shortbuf, link, 7);
	const int badpid=INT_MAX-1;//larger than pid_max
	assert_int_equal(get_job_progname(exe, PATH_MAX, badpid), 1);

	//memory usage increases after touching a large buffer
	size_t mem0=get_job_mem();
	assert_true(mem0>0);
	const size_t nbig=64<<20;
	char *big=(char *)malloc(nbig);
	volatile char *vbig=big;//prevent the compiler from removing the allocation
	for(size_t i=0; i<nbig; i+=1024) vbig[i]=1;
	size_t mem1=get_job_mem();
	assert_true(mem1>=mem0+(nbig>>10)/2);
	free(big);

	//launch time is after boot and before now
	double tlaunch=get_job_launchtime(getpid());
	assert_true(tlaunch>0);
#if defined(__linux__)
	double uptime=read_proc_num("/proc/uptime");
	assert_true(tlaunch<=uptime+1);
#endif
	assert_float_equal(get_job_launchtime(badpid), 0, 0);

	double load=get_usage_load();
	assert_true(load>=0);
#if defined(__linux__)
	assert_float_equal(load, read_proc_num("/proc/loadavg"), 5);
#endif
	double mem=get_usage_mem();
	assert_true(mem>0&&mem<1);

	long user1=-1, tot1=-1, user2=-1, tot2=-1;
	assert_int_equal(read_cpu_counter(&user1, &tot1), 0);
	assert_true(user1>=0&&tot1>0&&user1<=tot1);
	read_self_cpu();
	//keep the cpu busy for a while
	double t0=myclockd();
	volatile double sum=0;
	while(myclockd()<t0+0.3){
		for(int i=0; i<10000; i++) sum+=sqrt((double)i);
	}
	double self=read_self_cpu();
	assert_true(self>0.1&&self<=NCPU+0.5);
	assert_int_equal(read_cpu_counter(&user2, &tot2), 0);
	assert_true(user2>=user1&&tot2>tot1);

	double cpu=get_usage_cpu();
	assert_true(cpu>=0&&cpu<=1);
	assert_float_equal(get_usage_cpu(), cpu, 0);//cached within 0.1 seconds
	assert_true(get_usage_running()>=0);
	int avail=get_cpu_avail();
	assert_true(avail>=0&&avail<=NCPU);

	set_dirout("/");
	assert_string_equal(DIROUT, "/");
	set_dirout(SRCDIR);
	if(!mystrcmp(SRCDIR, HOME)){//HOME is replaced by ~
		char *expect2=stradd("~", &(SRCDIR)[strlen(HOME)], NULL);
		assert_string_equal(DIROUT, expect2);
		free(expect2);
	}else{
		assert_string_equal(DIROUT, SRCDIR);
	}
}
/** hashlittle() matches the lookup3 reference values for short data and samples long data. */
static void sysutil_hash(void **state){
	(void)state;
	uint32_t buf[64]={0};
	const char *str="Four score and seven years ago";
	memcpy(buf, str, strlen(str));//aligned copy
	//reference values from lookup3.c (driver5)
	assert_int_equal(hashlittle(buf, 0, 0), 0xdeadbeef);
	assert_int_equal(hashlittle(buf, 0, 0xdeadbeef), 0xbd5b7dde);
	assert_int_equal(hashlittle(buf, 30, 0), 0x17770551);
	assert_int_equal(hashlittle(buf, 30, 1), 0xcd628161);
	//bytes after length do not matter. Hash of every length is different.
	uint32_t hashes[25];
	uint32_t buf2[64];
	for(int len=0; len<25; len++){
		memcpy(buf2, buf, sizeof(buf));
		memset((char *)buf2+len, 0x5a, sizeof(buf2)-len);
		hashes[len]=hashlittle(buf, len, 3);
		assert_int_equal(hashlittle(buf2, len, 3), hashes[len]);
		for(int j=0; j<len; j++){
			assert_int_not_equal(hashes[j], hashes[len]);
		}
	}
	//long data: only a subset of words is used.
	const size_t nbyte=120000;//dstep=12: 3 words used out of every 36
	uint32_t *big=(uint32_t *)calloc(nbyte/4, sizeof(uint32_t));
	for(size_t i=0; i<nbyte/4; i++) big[i]=(uint32_t)(i*2654435761u);
	uint32_t h0=hashlittle(big, nbyte, 0);
	assert_int_equal(hashlittle(big, nbyte, 0), h0);
	assert_int_not_equal(hashlittle(big, nbyte, 1), h0);
	big[5]^=1;//not sampled
	assert_int_equal(hashlittle(big, nbyte, 0), h0);
	big[5]^=1;
	big[36]^=1;//sampled
	assert_int_not_equal(hashlittle(big, nbyte, 0), h0);
	free(big);
}
/** SHA1 and base64 encoding reproduce known test vectors. */
static void sysutil_sha1(void **state){
	(void)state;
	const struct{
		const char *in;
		const char *out;
	}b64[]={{"", ""}, {"f", "Zg=="}, {"fo", "Zm8="}, {"foo", "Zm9v"}, {"foob", "Zm9vYg=="}, {"fooba", "Zm9vYmE="}, {"foobar", "Zm9vYmFy"}};
	char out[64];
	for(unsigned int i=0; i<sizeof(b64)/sizeof(b64[0]); i++){
		base64_encode((const unsigned char *)b64[i].in, (int)strlen(b64[i].in), out);
		assert_string_equal(out, b64[i].out);
	}
	char a1000[1001];
	memset(a1000, 'a', 1000); a1000[1000]=0;
	char x60[61];
	memset(x60, 'x', 60); x60[60]=0;//padding needs an extra block
	const struct{
		const char *in;
		const char *out;
	}sha[]={
		{"", "2jmj7l5rSw0yVb/vlWAYkK/YBwk="},
		{"abc", "qZk+NkcGgWq6PiVxeFDCbJzQ2J0="},
		{"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq", "hJg+RBw70m66rkqh+VEp5eVGcPE="},
		{a1000, "KR6abGaZSUm1e6XmUDYemPw2sbo="},
		{x60, "Bs7S4HDljCxO2fK4y4kPDFEs5g0="},
		//websocket handshake example from RFC 6455
		{"dGhlIHNhbXBsZSBub25jZQ==258EAFA5-E914-47DA-95CA-C5AB0DC85B11", "s3pPLMBiTxaQ9kYGzzhZRbK+xOo="},
	};
	for(unsigned int i=0; i<sizeof(sha)/sizeof(sha[0]); i++){
		base64_sha1(sha[i].in, out);
		assert_string_equal(out, sha[i].out);
	}
	//incremental update gives the same digest. SHA1_CTX is opaque; use storage large enough for it.
	union{
		unsigned char c[512];
		uint64_t align;
	}ctxbuf;
	SHA1_CTX *ctx=(SHA1_CTX *)&ctxbuf;
	sha1_init(ctx);
	for(int i=0; i<10; i++){
		sha1_update(ctx, (const unsigned char *)a1000, 100);
	}
	unsigned char digest[20];
	sha1_final(ctx, digest);
	base64_encode(digest, 20, out);
	assert_string_equal(out, "KR6abGaZSUm1e6XmUDYemPw2sbo=");
	const unsigned char abc[20]={0xa9, 0x99, 0x3e, 0x36, 0x47, 0x06, 0x81, 0x6a, 0xba, 0x3e,
		0x25, 0x71, 0x78, 0x50, 0xc2, 0x6c, 0x9c, 0xd0, 0xd8, 0x9d};
	sha1_init(ctx);
	sha1_update(ctx, (const unsigned char *)"a", 1);
	sha1_update(ctx, (const unsigned char *)"bc", 2);
	sha1_final(ctx, digest);
	assert_memory_equal(digest, abc, 20);

	//session cookies are 64 hex characters and differ each time
	char *c1=generate_session_cookie();
	char *c2=generate_session_cookie();
	assert_non_null(c1);
	assert_non_null(c2);
	assert_int_equal(strlen(c1), 64);
	assert_int_equal(strspn(c1, "0123456789abcdef"), 64);
	assert_string_not_equal(c1, c2);
	free(c1);
	free(c2);
}
typedef struct{
	dlist list;//must be first
	int val;
}node_t;
/** Check that the list contains the values in order, in both directions. */
static void check_dlist(dlist *head, dlist *tail, const int *vals, int n){
	dlist *p=head;
	for(int i=0; i<n; i++){
		assert_non_null(p);
		assert_int_equal(((node_t *)p)->val, vals[i]);
		if(i==0) assert_null(p->prev);
		p=p->next;
	}
	assert_null(p);
	p=tail;
	for(int i=n-1; i>=0; i--){
		assert_non_null(p);
		assert_int_equal(((node_t *)p)->val, vals[i]);
		p=p->prev;
	}
	assert_null(p);
	if(!n){
		assert_null(head);
		assert_null(tail);
	}
}
/** The double linked list is a deque. */
static void sysutil_dlist(void **state){
	(void)state;
	node_t nodes[5];
	for(int i=0; i<5; i++){
		nodes[i].val=i;
		nodes[i].list.next=nodes[i].list.prev=(dlist *)&nodes[(i+1)%5];//garbage
	}
	dlist *head=NULL, *tail=NULL;
	assert_null(remove_head(&head, &tail));
	assert_null(remove_tail(&head, &tail));
	add_tail(&head, &tail, &nodes[1].list);
	add_tail(&head, &tail, &nodes[2].list);
	add_head(&head, &tail, &nodes[0].list);
	add_tail(&head, &tail, &nodes[3].list);
	{
		const int vals[]={0, 1, 2, 3};
		check_dlist(head, tail, vals, 4);
	}
	assert_ptr_equal(remove_head(&head, &tail), &nodes[0].list);
	assert_null(nodes[0].list.next);
	assert_ptr_equal(remove_tail(&head, &tail), &nodes[3].list);
	assert_null(nodes[3].list.prev);
	{
		const int vals[]={1, 2};
		check_dlist(head, tail, vals, 2);
	}
	add_head(&head, &tail, &nodes[4].list);
	{
		const int vals[]={4, 1, 2};
		check_dlist(head, tail, vals, 3);
	}
	assert_ptr_equal(remove_tail(&head, &tail), &nodes[2].list);
	assert_ptr_equal(remove_tail(&head, &tail), &nodes[1].list);
	assert_ptr_equal(remove_tail(&head, &tail), &nodes[4].list);
	check_dlist(head, tail, NULL, 0);
	add_head(&head, &tail, &nodes[2].list);//empty list
	assert_ptr_equal(head, tail);
	assert_ptr_equal(remove_head(&head, &tail), &nodes[2].list);
	check_dlist(head, tail, NULL, 0);
	//invalid arguments are ignored
	add_head(NULL, &tail, &nodes[0].list);
	add_head(&head, NULL, &nodes[0].list);
	add_tail(&head, &tail, NULL);
	add_tail(NULL, NULL, &nodes[0].list);
	check_dlist(head, tail, NULL, 0);
	assert_null(remove_head(NULL, &tail));
	assert_null(remove_tail(&head, NULL));
}

#if HAVE_LIBZSTD
typedef struct{
	int fd;
	const char *data;
	size_t len;
	int ans;
}zwriter_t;
/** Compress data into fd in pieces of varying size. */
static void *zstd_writer(void *arg){
	zwriter_t *w=(zwriter_t *)arg;
	zstd_t *zf=zstd_open(w->fd, 'w');
	w->ans=zf?0:-1;
	size_t pos=0;
	for(size_t step=1; zf&&pos<w->len; step=step*3+1){
		size_t n=pos+step>w->len?w->len-pos:step;
		w->ans|=zstd_write(zf, w->data+pos, n);
		pos+=n;
		if(step==13) w->ans|=zstd_flush(zf);
	}
	w->ans|=zstd_write(zf, w->data, 0);//nothing to write
	zstd_close(zf, 1);
	return NULL;
}
typedef struct{
	zstd_t *zw;
	int fd;
	const char *data;
	size_t len;
	int ans;
}ztrunc_t;
/** Write len bytes, flush, and shut down writing on fd. Run in a thread so the reader can drain the socket. */
static void *zstd_trunc_writer(void *arg){
	ztrunc_t *t=(ztrunc_t *)arg;
	t->ans=zstd_write(t->zw, t->data, t->len);
	t->ans|=zstd_flush(t->zw);
	shutdown(t->fd, SHUT_WR);
	return NULL;
}
#endif
/** zstd streams compress and decompress through a pipe. */
static void sysutil_zstd(void **state){
	(void)state;
#if HAVE_LIBZSTD
	alarm(30);
	const size_t len=1<<20;
	char *data=(char *)malloc(len);
	rand_t rstat;
	seed_rand(&rstat, 5);
	for(size_t i=0; i<len; i++){
		data[i]=(i&4096)?(char)(i%7):(char)lrand(&rstat);//mixed compressible and random data
	}
	int fd[2];
	assert_int_equal(pipe(fd), 0);
	zwriter_t w={fd[1], data, len, -100};
	pthread_t thread;
	pthread_create(&thread, NULL, zstd_writer, &w);
	zstd_t *zf=zstd_open(fd[0], 'r');
	assert_non_null(zf);
	char *data2=(char *)calloc(len, 1);
	size_t pos=0;
	for(size_t step=5; pos<len; step=step*2+3){
		size_t n=pos+step>len?len-pos:step;
		assert_int_equal(zstd_read(zf, data2+pos, n), 0);
		pos+=n;
	}
	assert_int_equal(zstd_read(zf, data2, 0), 0);
	char extra[4];
	assert_int_equal(zstd_read(zf, extra, 4), -2);//end of stream
	pthread_join(thread, NULL);
	assert_int_equal(w.ans, 0);
	assert_memory_equal(data, data2, len);
	//wrong direction
	assert_int_equal(zstd_write(zf, data, 4), -1);
	assert_int_equal(zstd_flush(zf), -1);
	assert_int_equal(zstd_close(zf, 0), 0);
	assert_int_not_equal(fcntl(fd[0], F_GETFD), -1);//still open
	close(fd[0]);

	//corrupted data
	assert_int_equal(pipe(fd), 0);
	assert_int_equal(stwrite(fd[1], data, 1000), 0);//random data
	close(fd[1]);
	zf=zstd_open(fd[0], 'r');
	assert_int_equal(zstd_read(zf, data2, 100), -1);
	zstd_close(zf, 1);

	//truncated data
	{
		int sv[2];
		assert_int_equal(socketpair(AF_UNIX, SOCK_STREAM, 0, sv), 0);
		zstd_t *zw=zstd_open(sv[0], 'w');
		assert_int_equal(zstd_read(zw, data2, 4), -1);//wrong direction
		//the socket buffer may be smaller than the compressed data (8 KiB in macOS), so write in a thread.
		ztrunc_t t={zw, sv[0], data, 100000, -100};
		pthread_t thread2;
		pthread_create(&thread2, NULL, zstd_trunc_writer, &t);
		zf=zstd_open(sv[1], 'r');
		assert_int_equal(zstd_read(zf, data2, 200000), -2);//only 100000 available
		pthread_join(thread2, NULL);
		assert_int_equal(t.ans, 0);
		assert_memory_equal(data, data2, 100000);
		zstd_close(zf, 1);
		//writing to a closed peer fails
		assert_int_equal(zstd_write(zw, data, len), -1);
		assert_int_equal(zstd_flush(zw), -1);
		zstd_close(zw, 1);
	}
	//reading a closed descriptor fails
	zf=zstd_open(fd[0], 'r');//already closed
	assert_int_equal(zstd_read(zf, data2, 4), -1);
	zstd_close(zf, 0);
	//invalid arguments
	assert_null(zstd_open(-1, 'r'));
	assert_null(zstd_open(0, 'x'));
	assert_int_equal(zstd_close(NULL, 0), -1);
	assert_int_equal(zstd_read(NULL, data2, 4), -1);
	assert_int_equal(zstd_write(NULL, data2, 4), -1);
	free(data);
	free(data2);
	alarm(0);
#else
	skip();
#endif
}

static unsigned int nmalloc=0, ncalloc=0, nrealloc=0, nfree=0;
static void *(*old_malloc)(size_t);
static void *(*old_calloc)(size_t, size_t);
static void *(*old_realloc)(void *, size_t);
static void (*old_free)(void *);
static void *count_malloc(size_t n){ atomic_add_fetch(&nmalloc, 1); return old_malloc(n); }
static void *count_calloc(size_t n, size_t m){ atomic_add_fetch(&ncalloc, 1); return old_calloc(n, m); }
static void *count_realloc(void *p, size_t n){ atomic_add_fetch(&nrealloc, 1); return old_realloc(p, n); }
static void count_free(void *p){ atomic_add_fetch(&nfree, 1); old_free(p); }
static void deinit_fun(void){}
/** Memory helpers: recalloc zeros new memory, verbose mode, external allocators, deinit registration, and environment. */
static void sysutil_mem(void **state){
	(void)state;
	//recalloc keeps old content and zeros the extension
	int *p=(int *)malloc_maos(4*sizeof(int));
	for(int i=0; i<4; i++) p[i]=i+1;
	p=(int *)recalloc_maos(p, 4*sizeof(int), 10*sizeof(int));
	for(int i=0; i<10; i++) assert_int_equal(p[i], i<4?i+1:0);
	p=(int *)recalloc_maos(p, 10*sizeof(int), 2*sizeof(int));//shrink
	assert_int_equal(p[1], 2);
	free_maos(p);

	//verbose mode prints each operation
	int verbose0=MEM_VERBOSE;
	MEM_VERBOSE=1;
	char *s=mycalloc(10, char);
	s=myrealloc(s, 20, char);
	char *s2=mymalloc(5, char);
	myfree(s2);
	myfree(s);
	MEM_VERBOSE=verbose0;

	//external allocators are used once registered
	old_malloc=malloc_default;
	old_calloc=calloc_default;
	old_realloc=realloc_default;
	old_free=free_default;
	register_malloc(NULL, NULL, NULL, NULL);//no effect
	assert_ptr_equal(malloc_default, old_malloc);
	assert_ptr_equal(free_default, old_free);
	register_malloc(count_malloc, count_calloc, count_realloc, count_free);
	void *m1=malloc_maos(16);
	void *m2=calloc_maos(2, 8);
	m1=realloc_maos(m1, 32);
	free_maos(m1);
	free_maos(m2);
	register_malloc(old_malloc, old_calloc, old_realloc, old_free);
	assert_true(nmalloc>=1);
	assert_true(ncalloc>=1);
	assert_true(nrealloc>=1);
	assert_true(nfree>=2);
	unsigned int nmalloc0=nmalloc;
	free_maos(malloc_maos(8));
	assert_int_equal(nmalloc, nmalloc0);

	//data registered with register_deinit() is freed at exit (verified by the memory accounting).
	register_deinit(deinit_fun, NULL);
	register_deinit(deinit_fun, NULL);//duplicate is ignored
	void *data=mymalloc(100, char);
	register_deinit(NULL, data);
	register_deinit(deinit_fun, data);//both are duplicates

	//environment variables
	int level0=LOG_LEVEL;
	setenv("MAOS_LOG_LEVEL", "-3", 1);
	read_sys_env();
	assert_int_equal(LOG_LEVEL, -3);
	unsetenv("MAOS_LOG_LEVEL");
	LOG_LEVEL=level0;
	read_sys_env();
	assert_int_equal(LOG_LEVEL, level0);//unchanged when not set
}
static void *thread_flag(void *arg){
	atomic_add_fetch((unsigned int *)arg, 1);
	return NULL;
}
static void *thread_mask(void *arg){
	thread_block_signal();
	sigset_t set;
	pthread_sigmask(SIG_BLOCK, NULL, &set);
	*(int *)arg=sigismember(&set, SIGINT)&&sigismember(&set, SIGTERM)&&sigismember(&set, SIGUSR1);
	return NULL;
}
static void *thread_fun_dummy(thread_t *data){
	return data;
}
/** Thread helpers create threads, block signals, and identify threads. */
static void sysutil_thread(void **state){
	(void)state;
	unsigned int flag=0;
	assert_true(thread_new(thread_flag, &flag)!=0);
	double t0=myclockd();
	while(!__atomic_load_n(&flag, __ATOMIC_SEQ_CST)&&myclockd()<t0+10){
		mysleep(0.001);
	}
	assert_int_equal(flag, 1);

	int blocked=0;
	pthread_t thread;
	pthread_create(&thread, NULL, thread_mask, &blocked);
	pthread_join(thread, NULL);
	assert_int_equal(blocked, 1);
	//the main thread is not affected
	sigset_t set;
	pthread_sigmask(SIG_BLOCK, NULL, &set);
	assert_false(sigismember(&set, SIGUSR1));
	block_signal(1);
	pthread_sigmask(SIG_BLOCK, NULL, &set);
	assert_true(sigismember(&set, SIGINT)&&sigismember(&set, SIGTERM));
	block_signal(0);
	pthread_sigmask(SIG_BLOCK, NULL, &set);
	assert_false(sigismember(&set, SIGINT)||sigismember(&set, SIGTERM));

#if _OPENMP
	assert_int_equal(thread_id(), 0);
	int nwrong=0;
#pragma omp parallel num_threads(4) reduction(+:nwrong)
	{
		if(thread_id()!=omp_get_thread_num()) nwrong++;
	}
	assert_int_equal(nwrong, 0);
#else
	assert_true(thread_id()>0);
#endif
	//function and data are set for active threads only
	int data=0;
	thread_t *thd=thread_prep(0, 10, 4, thread_fun_dummy, &data);
	assert_int_equal(thd[0].nthread, 4);//3 elements each, last one gets 1
	for(int i=0; i<4; i++){
		assert_ptr_equal(thd[i].fun, thread_fun_dummy);
		assert_ptr_equal(thd[i].data, &data);
		assert_ptr_equal(thd[i].fun(thd+i), thd+i);
	}
	assert_int_equal(thd[3].start, 9);
	assert_int_equal(thd[3].end, 10);
	free(thd);
	thd=thread_prep(0, 4, 8, thread_fun_dummy, &data);
	assert_int_equal(thd[0].nthread, 4);
	for(int i=0; i<8; i++){
		if(i<4){
			assert_ptr_equal(thd[i].fun, thread_fun_dummy);
			assert_int_equal(thd[i].end-thd[i].start, 1);
		}else{
			assert_null(thd[i].fun);
			assert_null(thd[i].data);
		}
	}
	free(thd);
}

UNITTEST_MAIN(
	cmocka_unit_test(sysutil_path),
	cmocka_unit_test(sysutil_process),
	cmocka_unit_test(sysutil_hash),
	cmocka_unit_test(sysutil_sha1),
	cmocka_unit_test(sysutil_dlist),
	cmocka_unit_test(sysutil_zstd),
	cmocka_unit_test(sysutil_mem),
	cmocka_unit_test(sysutil_thread),
)
