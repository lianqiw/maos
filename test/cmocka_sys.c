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
  Test utility routines in sys/. Replaces test_draw (print_backtrace) and
  test_plot (the shared memory interface it used no longer exists).
*/
#include <unistd.h>
#include <sys/stat.h>
#include "unittest.h"

/** print_backtrace() works without crashing. */
static void sys_backtrace(void **state){
	(void)state;
	print_backtrace();
}
/** mystrcmp() tests for a prefix and check_suffix() for a suffix. */
static void sys_prefix_suffix(void **state){
	(void)state;
	assert_int_equal(mystrcmp("abcdef", "abc"), 0);
	assert_int_equal(mystrcmp("abcdef", "abcdef"), 0);
	assert_int_not_equal(mystrcmp("abc", "abcdef"), 0);//a is shorter
	assert_int_not_equal(mystrcmp("abcdef", "abd"), 0);
	assert_int_not_equal(mystrcmp(NULL, "a"), 0);
	assert_int_not_equal(mystrcmp("a", NULL), 0);
	assert_int_not_equal(mystrcmp("", ""), 0);

	assert_int_equal(check_suffix("a.bin.gz", ".gz"), 1);
	assert_int_equal(check_suffix("a.bin.gz", ".bin.gz"), 1);
	assert_int_equal(check_suffix("a.bin.gz", ".bin"), 0);
	assert_int_equal(check_suffix("gz", ".gz"), 0);
	assert_int_equal(check_suffix(NULL, ".gz"), 0);
	assert_int_equal(check_suffix("a.gz", NULL), 0);
}
/** strscpy() always terminates and reports truncation; mystrndup() and mystrrep() handle lengths. */
static void sys_strcopy(void **state){
	(void)state;
	char buf[8];
	assert_int_equal(strscpy(buf, "abc", sizeof(buf)), 3);
	assert_string_equal(buf, "abc");
	assert_true(strscpy(buf, "abcdefghij", sizeof(buf))<0);//truncated
	assert_string_equal(buf, "abcdefg");
	assert_true(strscpy(buf, "abc", 0)<0);

	char *s=mystrndup("abcdef", 3);
	assert_string_equal(s, "abc");
	free(s);
	s=mystrndup("ab", 10);
	assert_string_equal(s, "ab");
	free(s);
	assert_null(mystrndup(NULL, 3));
	assert_null(mystrndup("abc", 0));

	char str[32]="/home/user/file";
	mystrrep(str, "/home/user", "~");
	assert_string_equal(str, "~/file");
	mystrrep(str, "/other", "~");//not a prefix: unchanged
	assert_string_equal(str, "~/file");
	mystrrep(str, "~", "/long");//substitute longer than prefix: unchanged
	assert_string_equal(str, "~/file");
}
/** mydirname(), stradd(), strnadd(), and argv2str() build the expected strings. */
static void sys_strbuild(void **state){
	(void)state;
	char *s=mydirname("/a/b/c.bin");
	assert_string_equal(s, "/a/b");
	free(s);
	s=mydirname("c.bin");
	assert_string_equal(s, ".");
	free(s);
	s=mydirname("/a/b/");//trailing / is ignored
	assert_string_equal(s, "/a");
	free(s);
	s=mydirname("/");//root
	assert_string_equal(s, "/");
	free(s);
	s=mydirname("/a");//file in root
	assert_string_equal(s, "/");
	free(s);
	s=mydirname("/a/");
	assert_string_equal(s, "/");
	free(s);
	assert_null(mydirname(NULL));
	assert_null(mydirname(""));

	s=stradd("ab", "cd", "", "ef", NULL);
	assert_string_equal(s, "abcdef");
	free(s);

	const char *argv[]={"x", "", NULL, "y"};
	s=strnadd(4, argv, ",");
	assert_string_equal(s, "x,y,");
	free(s);
	assert_null(strnadd(0, argv, ","));

	s=stradd(NULL, "ab", NULL);//first string may be NULL
	assert_string_equal(s, "ab");
	free(s);

	const char *argv2[]={"/usr/bin/maos", "-o", NULL, "out"};
	s=argv2str(4, argv2, " ");//prefixed by the starting directory
	char *expect=stradd(DIRSTART, "/maos -o out ", NULL);//keeps one trailing delimiter
	assert_string_equal(s, expect);
	free(expect);
	free(s);
	const char *argv3[]={"maos", "a", "b  "};
	s=argv2str(3, argv3, " ");//repeated trailing delimiters are trimmed to one
	expect=stradd(DIRSTART, "/maos a b ", NULL);
	assert_string_equal(s, expect);
	free(expect);
	free(s);
	assert_null(argv2str(0, argv2, " "));

	char buf[64];
	assert_int_equal(mysnprintf(buf, sizeof(buf), "%d-%s", 12, "ab"), 5);
	assert_string_equal(buf, "12-ab");
}
/** sec2str() picks the shortest representation of a duration. */
static void sys_sec2str(void **state){
	(void)state;
	char buf[32];
	sec2str(buf, sizeof(buf), 0);
	assert_string_equal(buf, "0");
	sec2str(buf, sizeof(buf), -5);
	assert_string_equal(buf, "0");
	sec2str(buf, sizeof(buf), 1.5);
	assert_string_equal(buf, "1.5s");
	sec2str(buf, sizeof(buf), 42);
	assert_string_equal(buf, "42s");
	sec2str(buf, sizeof(buf), 125);
	assert_string_equal(buf, "2m05s");
	sec2str(buf, sizeof(buf), 3*3600+7*60+30);
	assert_string_equal(buf, "3h07m");
	assert_int_equal(sec2str(buf, -1, 10), 0);
}
/** myclockd() is monotonic and mysleep() sleeps about the requested time. */
static void sys_clock(void **state){
	(void)state;
	double t0=myclockd();
	mysleep(0.05);
	double t1=myclockd();
	assert_true(t1-t0>=0.04);
	assert_true(t1-t0<5);
	assert_true(myclocki()>0);
}
/** thread_prep() splits a range into contiguous pieces that cover it exactly. */
static void sys_thread_prep(void **state){
	(void)state;
	int dummy=0;
	const long cases[][3]={{0, 100, 7}, {5, 9, 8}, {0, 1, 4}, {3, 3, 2}, {0, 64, 64}};
	for(unsigned int ic=0; ic<sizeof(cases)/sizeof(cases[0]); ic++){
		long start=cases[ic][0], end=cases[ic][1], nthread=cases[ic][2];
		thread_t *thd=thread_prep(start, end, nthread, NULL, &dummy);
		assert_non_null(thd);
		long nactive=thd[0].nthread;
		assert_true(nactive>=1&&nactive<=nthread);
		long pos=start;
		for(long it=0; it<nthread; it++){
			assert_int_equal(thd[it].ithread, it);
			if(it<nactive){
				assert_int_equal(thd[it].start, pos);
				assert_true(thd[it].end>=thd[it].start);
				pos=thd[it].end;
			}else{
				assert_int_equal(thd[it].end-thd[it].start, 0);
			}
		}
		assert_int_equal(pos, end);
		free(thd);
	}
	assert_null(thread_prep(0, 10, 0, NULL, NULL));
}
/** The atomic helpers are consistent when used from many threads. */
static void sys_atomic(void **state){
	(void)state;
	unsigned int count=0;
	unsigned int count2=100000;
#pragma omp parallel for
	for(int i=0; i<100000; i++){
		atomic_fetch_add(&count, 1);
		atomic_sub_fetch(&count2, 1);
	}
	assert_int_equal(count, 100000);
	assert_int_equal(count2, 0);
	assert_int_equal(atomic_add_fetch(&count, 5), 100005);
	assert_int_equal(atomic_fetch_sub(&count, 5), 100005);
	assert_int_equal(count, 100000);
}
#if defined(__linux__)
#define LINK_DIR "/proc/self" /*a symbolic link to a directory*/
#elif defined(__APPLE__)
#define LINK_DIR "/tmp" /*a symbolic link to private/tmp*/
#endif
/** File queries on existing paths. Operations that would modify files are only tried where they must fail. */
static void sys_files(void **state){
	(void)state;
	const char *fn=SRCDIR "/sys/misc.c";
	assert_int_equal(exist("/"), 1);
	assert_int_equal(exist(fn), 1);
#ifdef LINK_DIR
	assert_int_equal(exist(LINK_DIR), 1);
#endif
	assert_int_equal(exist("/nonexistent-maos-dir"), 0);
	assert_int_equal(exist(NULL), 0);
	assert_int_equal(isdir("/"), 1);
#ifdef LINK_DIR
	assert_int_equal(isdir(LINK_DIR), 1);//follows link
#endif
	assert_int_equal(isdir(fn), 0);
	assert_int_equal(isdir(NULL), 0);
	assert_int_equal(isfile(fn), 1);
	assert_int_equal(isfile("/"), 0);
	assert_int_equal(isfile(NULL), 0);
#ifdef LINK_DIR
	assert_int_equal(islink(LINK_DIR), 1);
#endif
	assert_int_equal(islink(fn), 0);
	assert_int_equal(islink(NULL), 0);
	struct stat buf;
	assert_int_equal(stat(fn, &buf), 0);
	assert_int_equal(flen(fn), buf.st_size);
	assert_true(flen(fn)>1000);
	assert_int_equal(flen("/nonexistent-maos-dir"), 0);
	assert_int_equal(flen(NULL), 0);
	assert_int_equal(fmtime(fn), buf.st_mtime);
	assert_int_equal(fmtime("/nonexistent-maos-dir"), 0);
	assert_int_equal(fmtime(NULL), 0);
	assert_true(available_space("/")>0);
	assert_int_equal(available_space("/nonexistent-maos-dir"), 0);
	//these fail without touching the file system
	assert_int_equal(copyfile("/nonexistent-maos-dir/a", "/nonexistent-maos-dir/b"), 1);
	assert_int_equal(copyfile(fn, "/nonexistent-maos-dir/b"), 1);
	assert_int_equal(mysymlink("/nonexistent-maos-dir/a", "/nonexistent-maos-dir/b"), -1);
	assert_int_not_equal(mylink("/nonexistent-maos-dir/a", "/nonexistent-maos-dir/b"), 0);
	assert_int_equal(exist("/nonexistent-maos-dir"), 0);
	mymkdir("%s/", "/dev");//already exists
	assert_int_equal(isdir("/dev"), 1);
	remove_file_older("/nonexistent-maos-dir", 1, 0);
	remove_file_older("/dev", 0, 1000000);//no .bin or .lock files
	rename_log(0, "nonexistent-maos-exe");//no log file to rename
}
/** print_file() writes the content of a file to stdout. */
static void sys_print_file(void **state){
	(void)state;
	const char *fn=SRCDIR "/AUTHORS";
	char expect[4096]={0};
	FILE *fp=fopen(fn, "r");
	assert_non_null(fp);
	size_t nexpect=fread(expect, 1, sizeof(expect)-1, fp);
	fclose(fp);
	assert_true(nexpect>0);
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	fflush(stdout);
	int saved=dup(1);
	dup2(fds[1], 1);
	print_file(fn);
	print_file("nonexistent-maos-file");//prints a warning only
	fflush(stdout);
	dup2(saved, 1);
	close(saved);
	close(fds[1]);
	char got[4096]={0};
	size_t ngot=0;
	ssize_t n;
	while(ngot<sizeof(got)-1&&(n=read(fds[0], got+ngot, sizeof(got)-1-ngot))>0) ngot+=(size_t)n;
	close(fds[0]);
	assert_int_equal(ngot, nexpect);
	assert_memory_equal(got, expect, nexpect);
}
/** myabspath() resolves ., .., and ~; expand_filename() expands ~. */
static void sys_abspath(void **state){
	(void)state;
	char *cwd=getcwd(NULL, 0);
	char *s=myabspath(NULL);
	assert_string_equal(s, cwd);
	free(s);
	s=myabspath(".");
	assert_string_equal(s, cwd);
	free(s);
	s=myabspath("/a/b");
	assert_string_equal(s, "/a/b");
	free(s);
	s=myabspath("./abc");
	char *expect=stradd(cwd, "/abc", NULL);
	assert_string_equal(s, expect);
	free(s);
	s=myabspath("abc");
	assert_string_equal(s, expect);
	free(s);
	free(expect);
	s=myabspath("~/abc");
	expect=stradd(HOME, "/abc", NULL);
	assert_string_equal(s, expect);
	free(s);
	free(expect);
	s=myabspath("../");
	char *parent=mydirname(cwd);
	assert_string_equal(s, parent);
	free(s);
	s=myabspath("..//abc");
	expect=stradd(parent, "/abc", NULL);
	assert_string_equal(s, expect);
	free(s);
	free(expect);
	s=myabspath("../../x");
	char *parent2=mydirname(parent);
	expect=stradd(parent2, "/x", NULL);
	assert_string_equal(s, expect);
	free(s);
	free(expect);
	free(parent2);
	free(parent);
	free(cwd);

	s=expand_filename("~/abc");
	expect=stradd(HOME, "/abc", NULL);
	assert_string_equal(s, expect);
	free(s);
	free(expect);
	s=expand_filename("/abc");
	assert_string_equal(s, "/abc");
	free(s);
}
/** myasctime() and strtime_pid() format the time. */
static void sys_time(void **state){
	(void)state;
	time_t t=1700000000;
	struct tm tm;
	localtime_r(&t, &tm);
	char expect[64];
	snprintf(expect, sizeof(expect), "%04d/%02d/%02d %02d:%02d:%02d", tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	assert_string_equal(myasctime(t), expect);
	assert_int_equal(strlen(myasctime(0)), 19);
	char *s=strtime_pid();
	char pid[32];
	snprintf(pid, sizeof(pid), "-%05d", (int)getpid());
	assert_true(check_suffix(s, pid));
	assert_int_equal(strlen(s), strlen("2023-11-14-221320")+strlen(pid));
	int year=0, mon=0, day=0;
	assert_int_equal(sscanf(s, "%d-%d-%d", &year, &mon, &day), 3);
	assert_true(year>=2023&&mon>=1&&mon<=12&&day>=1&&day<=31);
	free(s);
}
/** Miscellaneous string helpers. */
static void sys_strmisc(void **state){
	(void)state;
	assert_null(mystrdup(NULL));
	char *s=mystrdup("abc");
	assert_string_equal(s, "abc");
	free(s);
	char buf[4];
	const char *longstr="abcdef";
	assert_int_equal(mysnprintf(buf, sizeof(buf), "%s", longstr), 6);//truncated
	assert_string_equal(buf, "abc");
	free_strarr(NULL, 3);
	char **arr=mycalloc(2, char *);
	arr[0]=strdup("a");
	free_strarr(arr, 2);//NULL entry is fine
	assert_int_equal(default_color(0), 0x0000FF);
	assert_int_equal(default_color(12), default_color_table[1]);
	print_version();
	//argv2str with default delimiter replaces newline in arguments
	const char *argv[]={"maos", "a\nb", "c"};
	s=argv2str(3, argv, NULL);
	char *expect=stradd(DIRSTART, "/maos a b c ", NULL);
	assert_string_equal(s, expect);
	free(expect);
	free(s);
	//newline as delimiter is kept
	s=argv2str(3, argv, "\n");
	expect=stradd(DIRSTART, "/maos\na\nb\nc\n", NULL);
	assert_string_equal(s, expect);
	free(expect);
	free(s);
	//repeated trailing delimiters are trimmed to one
	const char *argv3[]={"a", "b,"};
	s=strnadd(2, argv3, ",");
	assert_string_equal(s, "a,b,");
	free(s);
	char str[16]="abcdef";
	mystrrep(str, NULL, "x");
	mystrrep(NULL, "a", "x");
	mystrrep(str, "abcdefgh", "x");//prefix longer than string
	assert_string_equal(str, "abcdef");
	mystrrep(str, "abc", "XY");
	assert_string_equal(str, "XYdef");
	char tmp[32];
	sec2str(tmp, sizeof(tmp), 0.25);
	assert_string_equal(tmp, "0.25s");
	sec2str(tmp, sizeof(tmp), 10);
	assert_string_equal(tmp, "10s");
	sec2str(tmp, sizeof(tmp), 3600);
	assert_string_equal(tmp, "1h00m");
}
/** mypause() returns according to the key read from the descriptor. */
static void sys_mypause(void **state){
	(void)state;
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	const char keys[]="c7x";
	const int expect[]={0, 7, 1};
	for(int i=0; i<3; i++){
		assert_int_equal(write(fds[1], keys+i, 1), 1);
		assert_int_equal(mypause(-1, fds[0]), expect[i]);
	}
	assert_int_equal(write(fds[1], "5", 1), 1);
	assert_int_equal(mypause(fds[0], fds[0]), 5);
	close(fds[1]);
	assert_int_equal(mypause(fds[0], -1), 1);//end of file
	close(fds[0]);
}
void str_swap(char **pdest, char *begin, char *end);
/** str_swap() moves a segment to the destination or removes it. */
static void sys_str_swap(void **state){
	(void)state;
	char str[32]="0123456789";
	char *dest=str;
	str_swap(&dest, str+4, str+7);//move 456 to the front
	assert_string_equal(str, "4560123789");
	assert_ptr_equal(dest, str+3);
	str_swap(&dest, str+3, str+5);//at the destination: just advance
	assert_string_equal(str, "4560123789");
	assert_ptr_equal(dest, str+5);
	char *end=str+strlen(str);
	str_swap(&end, str+2, str+4);//remove 60
	assert_string_equal(str, "45123789");
	assert_ptr_equal(end, str+8);
}
static int fun_int_val=0;
static double fun_dbl_val=0;
static char fun_str_val[32]={0};
static void fun_int(int val){ fun_int_val=val; }
static void fun_dbl(double val){ fun_dbl_val=val; }
static void fun_str(const char *val){ strscpy(fun_str_val, val, sizeof(fun_str_val)); }
/** parse_argopt() extracts options and leaves configuration entries one per line. */
static void sys_parse_argopt(void **state){
	(void)state;
	int detach=0, nthread=0, ngpu=0, ndarr=0, flag=0;
	double scale=0;
	char *dir=NULL, *name=NULL;
	int *gpus=NULL;
	double *darr=NULL;
	argopt_t options[]={
		{"help", 'h', M_INT, 0, 1, (void *)fun_int, NULL},
		{"detach", 'd', M_INT, 0, 0, &detach, NULL},
		{"nthread", 'n', M_INT, 1, 0, &nthread, NULL},
		{"output", 'o', M_STR, 1, 0, &dir, NULL},
		{"gpu", 'g', M_INT, 2, 0, &gpus, &ngpu},
		{"scale", 's', M_DBL, 1, 0, &scale, NULL},
		{"arr", 'a', M_DBL, 2, 0, &darr, &ndarr},
		{"fdbl", 'f', M_DBL, 1, 1, (void *)fun_dbl, NULL},
		{"fstr", 'F', M_STR, 1, 1, (void *)fun_str, NULL},
		{"name", 'N', M_STR, 1, 0, &name, NULL},
		{"flag", 'l', M_INT, 1, 0, &flag, NULL},
		{"skip", 'x', 0, 0, 0, NULL, NULL},
		{NULL, 0, 0, 0, 0, NULL, NULL}
	};
	char *cmds=strdup("-d -n 4 --output=out1 -g1 -g 2 --scale 0.5 -a 1.5 -a2 sim.conf"
		" key=value -h --fdbl=3 -F 'quoted str' -x -z extra -o out2 dm.dx = [1 2 3] str='a b' --name");
	parse_argopt(cmds, options);
	assert_int_equal(detach, 1);
	assert_int_equal(nthread, 4);
	assert_string_equal(dir, "out2");
	assert_int_equal(ngpu, 2);
	assert_int_equal(gpus[0], 1);
	assert_int_equal(gpus[1], 2);
	assert_float_equal(scale, 0.5, 0);
	assert_int_equal(ndarr, 2);
	assert_float_equal(darr[0], 1.5, 0);
	assert_float_equal(darr[1], 2, 0);
	assert_int_equal(fun_int_val, 1);
	assert_float_equal(fun_dbl_val, 3, 0);
	assert_string_equal(fun_str_val, "quoted str");
	assert_string_equal(name, "Unknown");//no value given
	assert_int_equal(flag, 0);
	//options are removed. Unknown options, conf files and key=value entries remain.
	assert_string_equal(cmds, "      \nsim.conf\nkey=value     -z extra\ndm.dx = [1 2 3]\nstr='a b' ");
	assert_non_null(strstr(cmds, "sim.conf\n"));
	assert_non_null(strstr(cmds, "\nkey=value"));
	assert_non_null(strstr(cmds, "-z extra"));
	assert_non_null(strstr(cmds, "\ndm.dx = [1 2 3]"));
	assert_non_null(strstr(cmds, "\nstr='a b'"));
	assert_null(strstr(cmds, "-d"));
	assert_null(strstr(cmds, "out1"));
	assert_null(strstr(cmds, "quoted"));
	assert_null(strstr(cmds, "-x"));
	free(cmds);
	free(dir);
	free(name);
	free(gpus);
	free(darr);
	//NULL is ignored. Without options, only the configuration is split.
	parse_argopt(NULL, options);
	cmds=strdup("a.conf b=1 -o dir 'c' \"d=2\"");
	parse_argopt(cmds, NULL);
	assert_string_equal(cmds, "a.conf\nb=1\n-o dir 'c' \"d=2\"");
	free(cmds);
	//nested brackets are kept together. An unclosed quote is tolerated.
	cmds=strdup("x=[[1 2] 3] y='abc");
	parse_argopt(cmds, NULL);
	assert_string_equal(cmds, "x=[[1 2] 3]\ny='abc");
	free(cmds);
	//invalid usage reports error and continues safely.
	argopt_t options2[]={
		{"gpu", 'g', M_INT, 2, 0, &gpus, &ngpu},
		{"arr", 'a', M_DBL, 2, 0, &darr, &ndarr},
		{"farr", 'f', M_INT, 2, 1, (void *)fun_int, NULL},//array to function is not supported
		{"name", 'N', M_STR, 1, 0, &name, NULL},
		{NULL, 0, 0, 0, 0, NULL, NULL}
	};
	ngpu=0; ndarr=0; gpus=NULL; darr=NULL; name=NULL;
	cmds=strdup("-f 3 -N 'unterminated -g");//-g without value
	parse_argopt(cmds, options2);
	assert_int_equal(ngpu, 0);
	assert_null(gpus);
	assert_string_equal(name, "unterminated -g");//takes the remaining string
	free(name); name=NULL;
	free(cmds);
	cmds=strdup("-a");//-a without value
	parse_argopt(cmds, options2);
	assert_int_equal(ndarr, 0);
	assert_null(darr);
	free(cmds);
}

UNITTEST_MAIN(
	cmocka_unit_test(sys_backtrace),
	cmocka_unit_test(sys_prefix_suffix),
	cmocka_unit_test(sys_strcopy),
	cmocka_unit_test(sys_strbuild),
	cmocka_unit_test(sys_sec2str),
	cmocka_unit_test(sys_clock),
	cmocka_unit_test(sys_thread_prep),
	cmocka_unit_test(sys_atomic),
	cmocka_unit_test(sys_files),
	cmocka_unit_test(sys_print_file),
	cmocka_unit_test(sys_abspath),
	cmocka_unit_test(sys_time),
	cmocka_unit_test(sys_strmisc),
	cmocka_unit_test(sys_mypause),
	cmocka_unit_test(sys_str_swap),
	cmocka_unit_test(sys_parse_argopt),
)
