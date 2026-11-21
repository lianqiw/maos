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
  Tests of the configuration parser lib/readcfg.c. Entries are supplied as
  key=value strings, the same way maos passes command line overrides. A few
  tests read existing .conf files from the source tree (read only). The
  parser has global state, so each test ends with close_config(), which also
  checks for unused entries.
*/
#include "unittest.h"

#define CONFDIR SRCDIR "/config/maos"

/*close the configuration without saving and return whether error() was called.*/
static int close_check(void){
	signal_caught=0;
	close_config("%s", "");
	int ans=signal_caught==SIGUSR2;
	signal_caught=0;
	return ans;
}

/** Scalars, strings, arrays and matrices supplied in a single string. */
static void cfg_basic(void **state){
	(void)state;
	open_config("sim.a=3 sim.b=2.5e-3 sim.s='hello world' sim.e= powfs1.n=7 "
		"arr=[1 2 3] darr=[1.5,2.5] mat=[1 2 3;4 5 6] sarr=[\"ab\",\"c d\",e] sarr2=['x y'] expr=2*[1 2]+1 neg=-4", 0);
	assert_int_equal(readcfg_peek("sim.a"), 1);
	assert_int_equal(readcfg_peek("sim.none"), 0);
	assert_int_equal(readcfg_peek_n("arr"), 3);
	assert_int_equal(readcfg_peek_n("darr"), 2);
	assert_int_equal(readcfg_peek_n("sim.a"), 1);
	//readcfg_peek_n() only recognizes double quoted string arrays.
	assert_int_equal(readcfg_peek_n("sarr"), 3);
	assert_int_equal(readcfg_peek_n("sim.none"), 0);
	//direct string entries have priority 1
	assert_int_equal(readcfg_peek_priority("sim.a"), 1);
	assert_int_equal(readcfg_peek_priority("sim.none"), 0);
	assert_true(readcfg_peek_index("sim.b")>readcfg_peek_index("sim.a"));
	assert_int_equal(readcfg_peek_index("sim.none"), 0);

	assert_int_equal(readcfg_int("sim.a"), 3);
	assert_float_equal(readcfg_dbl("sim.b"), 2.5e-3, 1e-15);
	assert_int_equal(readcfg_int("powfs%d.n", 1), 7);
	assert_int_equal(readcfg_int("neg"), -4);
	char *str=readcfg_str("sim.s");
	assert_string_equal(str, "hello world");
	free(str);
	assert_null(readcfg_str("sim.e"));//empty value
	int *iarr=NULL;
	assert_int_equal(readcfg_intarr(&iarr, 0, 0, "arr"), 3);
	assert_int_equal(iarr[0], 1);
	assert_int_equal(iarr[2], 3);
	free(iarr);
	real *darr=NULL;
	assert_int_equal(readcfg_dblarr(&darr, 2, 0, "darr"), 2);
	assert_float_equal(darr[1], 2.5, 0);
	free(darr);
	//columns are separated by semi-colon
	dmat *m=readcfg_dmat(0, 0, "mat");
	assert_int_equal(NX(m), 3);
	assert_int_equal(NY(m), 2);
	assert_float_equal(P(m, 0, 1), 4, 0);
	assert_float_equal(P(m, 2, 0), 3, 0);
	dfree(m);
	lmat *l=readcfg_lmat(0, 0, "expr");
	assert_int_equal(PN(l), 2);
	assert_int_equal(P(l, 0), 3);
	assert_int_equal(P(l, 1), 5);
	lfree(l);
	char **sarr=NULL;
	int ns=readcfg_strarr(&sarr, 0, 0, "sarr");
	assert_int_equal(ns, 3);
	assert_string_equal(sarr[0], "ab");
	assert_string_equal(sarr[1], "c d");
	assert_string_equal(sarr[2], "e");
	for(int i=0; i<ns; i++) free(sarr[i]);
	free(sarr);
	//preallocated with relax
	char *sarr2[2]={NULL, NULL};
	char **psarr2=sarr2;
	assert_int_equal(readcfg_strarr(&psarr2, 2, 1, "sarr2"), 2);
	assert_string_equal(sarr2[0], "x y");
	assert_string_equal(sarr2[1], "x y");
	free(sarr2[0]); free(sarr2[1]);
	//missing keys
	assert_null(readcfg_dmat(0, 0, "sim.none"));
	assert_null(readcfg_lmat(0, 0, "sim.none"));
	assert_null(readcfg_str("sim.none"));
	assert_int_equal(readcfg_int("sim.none"), 0);
	assert_float_equal(readcfg_dbl("sim.none"), 0, 0);
	readcfg_ignore("sim.none");//no error
	//all keys are used
	assert_int_equal(close_check(), 0);
	//after close, nothing is left
	assert_int_equal(readcfg_peek("sim.a"), 0);
}

/** relax: fewer entries are filled with the last value. */
static void cfg_relax(void **state){
	(void)state;
	open_config("one=[2] two=[1 2] one2=[5] two2=[3 4] one3=7 two3=[1 2] two4=[1 2]", 0);
	real *d=NULL;
	assert_int_equal(readcfg_dblarr(&d, 3, 1, "one"), 3);
	assert_float_equal(d[0], 2, 0);
	assert_float_equal(d[2], 2, 0);
	free(d);
	int *i=NULL;
	assert_int_equal(readcfg_intarr(&i, 4, 2, "two"), 4);
	assert_int_equal(i[1], 2);
	assert_int_equal(i[3], 2);
	free(i);
	dmat *m=readcfg_dmat(3, 1, "one2");
	assert_int_equal(PN(m), 3);
	assert_float_equal(P(m, 2), 5, 0);
	dfree(m);
	m=readcfg_dmat(4, 2, "two2");
	assert_int_equal(PN(m), 4);
	assert_float_equal(P(m, 3), 4, 0);
	dfree(m);
	lmat *l=readcfg_lmat(2, 1, "one3");
	assert_int_equal(PN(l), 2);
	assert_int_equal(P(l, 1), 7);
	lfree(l);
	//relax=1 does not allow 2 of 3 entries.
	signal_caught=0;
	d=NULL;
	readcfg_dblarr(&d, 3, 1, "two3");
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	free(d);
	//relax=0 needs all entries
	l=readcfg_lmat(3, 0, "two4");
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	lfree(l);
	assert_int_equal(close_check(), 0);
}

/** Overrides follow priority and input order. */
static void cfg_priority(void **state){
	(void)state;
	open_config("x=1 y=[1]", 0);//priority 1
	open_config("x=2", 0);//same priority, later: overrides
	open_config("x=3", -1);//lower priority: ignored
	assert_int_equal(readcfg_peek_priority("x"), 1);
	open_config("y=[5]", 4);//higher priority
	assert_int_equal(readcfg_peek_priority("y"), 5);
	open_config("y=[6]", 0);//lower priority: ignored
	open_config("y=[5]", 0);//same value: nothing changes
	assert_int_equal(readcfg_peek_priority("y"), 5);
	assert_int_equal(readcfg_int("x"), 2);
	dmat *y=readcfg_dmat(1, 0, "y");
	assert_float_equal(P(y, 0), 5, 0);
	dfree(y);
	//leading slashes and spaces are skipped. Empty strings are ignored.
	open_config(" //z=9", 0);
	open_config("   ", 0);
	open_config("", 0);
	open_config(NULL, 0);
	assert_int_equal(readcfg_int("z"), 9);
	assert_int_equal(close_check(), 0);
}

/** += appends to and -= removes from arrays. */
static void cfg_append(void **state){
	(void)state;
	open_config("v=[1 2] w=[] u=[1 2 3 4] s=abc t=[1]", 0);
	open_config("v+=[3 4]", 0);
	open_config("w+=[5]", 0);
	open_config("u-=[3 4]", 0);
	dmat *v=readcfg_dmat(0, 0, "v");
	assert_int_equal(PN(v), 4);
	assert_float_equal(P(v, 3), 4, 0);
	dfree(v);
	dmat *w=readcfg_dmat(0, 0, "w");
	assert_int_equal(PN(w), 1);
	assert_float_equal(P(w, 0), 5, 0);
	dfree(w);
	dmat *u=readcfg_dmat(0, 0, "u");
	assert_int_equal(PN(u), 2);
	assert_float_equal(P(u, 1), 2, 0);
	dfree(u);
	//invalid append or removal
	signal_caught=0;
	open_config("s+=[1]", 0);//old value is not an array
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	open_config("t+=2", 0);//new value is not an array
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	open_config("t-=[9]", 0);//not at the end
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	char *s=readcfg_str("s");
	assert_string_equal(s, "abc");
	free(s);
	dmat *t=readcfg_dmat(0, 0, "t");
	assert_int_equal(PN(t), 1);
	dfree(t);
	assert_int_equal(close_check(), 0);
}

/** Renamed and equivalent keys, removed keys, and references. */
static void cfg_rename(void **state){
	(void)state;
	//change records as in changes.conf
	open_config("old.key=>new.key alt.key==main.key gone.key==(removed)", -1);
	open_config("old.key=7 alt.key=[1 2] gone.key=3 other=ref:main.key", 0);
	assert_int_equal(readcfg_peek("gone.key"), 1);//the change record itself
	assert_int_equal(readcfg_int("new.key"), 7);
	assert_int_equal(readcfg_peek_n("main.key"), 2);
	//reference to another key
	dmat *o=readcfg_dmat(0, 0, "other");
	assert_int_equal(PN(o), 2);
	assert_float_equal(P(o, 1), 2, 0);
	dfree(o);
	dmat *mk=readcfg_dmat(0, 0, "main.key");
	assert_int_equal(PN(mk), 2);
	dfree(mk);
	//using a change record is an error
	signal_caught=0;
	readcfg_int("old.key");
	assert_int_equal(signal_caught, SIGUSR2);
	assert_int_equal(close_check(), 0);
	//a renamed key given an empty value
	open_config("old.key=>new.key", -1);
	open_config("old.key=", 0);
	assert_int_equal(readcfg_peek("new.key"), 1);
	assert_null(readcfg_str("new.key"));
	assert_int_equal(close_check(), 0);
}

/** Unknown, unused or multiply used entries are reported. */
static void cfg_errors(void **state){
	(void)state;
	//unused entry with priority 0 (default) is only a warning
	open_config("unused.key=1", -1);
	assert_int_equal(close_check(), 0);
	//dbg. entries are never an error
	open_config("dbg.something=1", 0);
	assert_int_equal(close_check(), 0);
	//unused entry from the command line is an error
	open_config("unused.key=1", 0);
	assert_int_equal(close_check(), 1);
	//used twice
	open_config("a=1", 0);
	assert_int_equal(readcfg_int("a"), 1);
	signal_caught=0;
	readcfg_int("a");
	assert_int_equal(signal_caught, SIGUSR2);
	assert_int_equal(close_check(), 0);//already reported
	//invalid values
	open_config("i=1.5 d=abc", 0);
	signal_caught=0;
	readcfg_int("i");
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	readcfg_dbl("d");
	assert_int_equal(signal_caught, SIGUSR2);
	assert_int_equal(close_check(), 0);
	//invalid lines
	signal_caught=0;
	open_config("garbage", 0);
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	open_config("=5", 0);
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	//lines starting with __ are skipped
	open_config("__reset__", 0);
	assert_int_equal(signal_caught, 0);
	//missing file
	open_config("nonexistent_unittest.conf", -1);//silent
	assert_int_equal(signal_caught, 0);
	open_config("nonexistent_unittest.conf", 0);
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	assert_int_equal(close_check(), 0);
}

/** Prefixed entries take the priority and index of an existing key. */
static void cfg_prefix(void **state){
	(void)state;
	open_config("powfs.order=[4]", 3);
	open_config_prefix("order=8 extra=2", "powfs0_", "powfs.order");
	//string entries get one higher priority than given
	assert_int_equal(readcfg_peek_priority("powfs0_order"), readcfg_peek_priority("powfs.order")+1);
	assert_int_equal(readcfg_peek_index("powfs0_order"), readcfg_peek_index("powfs.order"));
	assert_int_equal(readcfg_int("powfs0_order"), 8);
	assert_int_equal(readcfg_int("powfs0_extra"), 2);
	lmat *order=readcfg_lmat(1, 0, "powfs.order");
	assert_int_equal(P(order, 0), 4);
	lfree(order);
	//a missing file with prefix is only a warning
	signal_caught=0;
	open_config_prefix("nonexistent_unittest.conf", "powfs1_", "powfs.order");
	assert_int_equal(signal_caught, 0);
	assert_int_equal(close_check(), 0);
}

/** Read existing .conf files (read only): include, path, changes.conf, and environment variables. */
static void cfg_files(void **state){
	(void)state;
	//absolute file name with include= lines
	open_config(CONFDIR "/atmr_single.conf", 0);
	dmat *ht=readcfg_dmat(0, 0, "atmr.ht");
	assert_int_equal(PN(ht), 1);
	assert_float_equal(P(ht, 0), 0, 0);
	dfree(ht);
	//the file has priority 0
	assert_int_equal(readcfg_peek_priority("atmr.wt"), 0);
	readcfg_ignore("atmr.wt");
	//the directory of the default config is added to the search path, so include= works.
	//(A string that ends with .conf is treated as a file name.)
	//file entries can be overriden by the command line
	open_config("include=dm_single.conf dm.offset=[0.5]", 0);
	assert_int_equal(readcfg_peek_n("dm.ht"), 1);
	assert_int_equal(readcfg_peek_priority("dm.ht"), 0);
	dmat *off=readcfg_dmat(1, 0, "dm.offset");
	assert_float_equal(P(off, 0), 0.5, 0);
	dfree(off);
	readcfg_ignore("dm.ht");
	assert_int_equal(close_check(), 0);

	//path= adds a search directory. changes.conf is then found and loaded on the first call.
	open_config("path=" CONFDIR, 0);
	assert_int_equal(close_check(), 0);
	open_config("atm.zadeg=30", 0);//renamed to sim.zadeg in changes.conf
	assert_int_equal(readcfg_int("sim.zadeg"), 30);
	open_config("fov_sq34.conf", 0);//nested includes
	assert_int_equal(readcfg_peek_n("evl.thetax"), 9);
	assert_int_equal(readcfg_peek_priority("evl.thetax"), 0);
	dmat *wt=readcfg_dmat(9, 0, "evl.wt");
	assert_float_equal(dsum(wt), 1, 1e-12);
	dfree(wt);
	//MAOS_ entries set environment variables
	open_config("MAOS_UNITTEST_DUMMY=5", 0);
	LOG_LEVEL=-4;//read_sys_env() may reset it
	assert_non_null(getenv("MAOS_UNITTEST_DUMMY"));
	assert_string_equal(getenv("MAOS_UNITTEST_DUMMY"), "5");
	unsetenv("MAOS_UNITTEST_DUMMY");
	assert_int_equal(close_check(), 0);//entries from files have priority 0: unused entries are warnings
	rmpath(CONFDIR);
}

/*Clear the configuration even if a test fails.*/
static int teardown(void **state){
	(void)state;
	close_config("%s", "");
	signal_caught=0;
	return 0;
}
#define TEST(f) cmocka_unit_test_teardown(f, teardown)
UNITTEST_MAIN(
	TEST(cfg_basic),
	TEST(cfg_relax),
	TEST(cfg_priority),
	TEST(cfg_append),
	TEST(cfg_rename),
	TEST(cfg_errors),
	TEST(cfg_prefix),
	TEST(cfg_files)
)
