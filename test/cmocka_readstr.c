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
#include "unittest.h"

static void readstr_basic(void** state){
	(void)state;
	assert_float_equal(readstr_num("test","1+2*3", NULL), 7, 0);
	assert_float_equal(readstr_num("test","1*2*3", NULL), 6, 0);
	assert_float_equal(readstr_num("test","1*2+3", NULL), 5, 0);
	assert_float_equal(readstr_num("test","1+2+3", NULL), 6, 0);
	assert_float_equal(readstr_num("test","1+2-3", NULL), 0, 0);
	assert_float_equal(readstr_num("test","1+2/3", NULL), 1+2./3, 0.01);
	assert_float_equal(readstr_num("test","1+2+3+4*5", NULL), 1+2+3+4*5, 0.01);
	assert_float_equal(readstr_num("test", "1,", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1;", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1 1", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1;1", NULL), 1, 0.01);
	assert_true(isnan(readstr_num("test","/0.5+2+3+4*5", NULL)));
	assert_true(isnan(readstr_num("test","*0.5+2+3+4*5", NULL)));
	assert_true(!isnan(readstr_num("test","+0.5+2+3+4*5", NULL)));
}
static void readstr_array(void **state){
	(void)state;
	double *ret=0;
	int nrow;
	int ncol;

	//Simple read
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "[]"), 0);
	assert_int_equal(nrow, 0);
	assert_int_equal(ncol, 0);
	free(ret);ret=NULL;

	//Simple read
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2]"), 2);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1],1,0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 1);
	free(ret);ret=NULL;

	//Simple read of a matrix
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2; 3 4; 5 6]"), 6);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1], 1, 0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 3);
	free(ret);ret=NULL;

	//Simple read of a matrix; invalid entry
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2; 3 4; 5 ]"), -1);
	free(ret);ret=NULL;

	//Simple read without bracket. Ignore after ;
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2 1;3"), 2);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1], 1, 0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 1);
	free(ret);ret=NULL;

	//Relaxed read type 1
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 10, 1, M_DBL, "test", "2/[2]"), 10);
	assert_float_equal(ret[1],1,0.0001);
	assert_float_equal(ret[9],1,0.0001);
	assert_int_equal(nrow,10);
	assert_int_equal(ncol,1);
	free(ret);ret=NULL;
	//Relaxed read type 2
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 10, 2, M_DBL, "test", "2/[1 2]"), 10);
	assert_float_equal(ret[2],1,0.0001);
	assert_float_equal(ret[9],1,0.0001);
	assert_int_equal(nrow,10);
	assert_int_equal(ncol,1);
	free(ret);ret=NULL;

	//Marix
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "3.2/[1/2*3+2*4/2 1*2 1+1*2 1-2; 1 2 3 -4]*2+1*2"), 8);
	assert_int_equal(nrow, 4);
	assert_int_equal(ncol, 2);
	assert_float_equal(ret[0], 3.2/(1./2*3+2.*4/2)*2+1.*2, 0.0001);
	assert_float_equal(ret[7], 3.2/(-4)*2+1.*2, 0.0001);
	free(ret);ret=NULL;

	//If not enough data is available, a relaxed vector is read if relax is set but not 1. matrix dimension is not honored.
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 9, 2, M_DBL, "test", "3.2/[1/2*3+2*4/2 1*2 1+1*2 1-2; 1 2 3 -4]*2+1*2"), 9);
	assert_int_equal(nrow,9);
	assert_int_equal(ncol,1);
	assert_float_equal(ret[8],ret[7],EPS);
	free(ret);ret=NULL;
}

/** readstr_num() handles operator precedence, end pointer, and invalid input. */
static void readstr_num_edge(void **state){
	(void)state;
	char *endptr=NULL;
	const char *str="  1.5e2*2, next";
	assert_float_equal(readstr_num("test", str, &endptr), 300, 0);
	assert_string_equal(endptr, ", next");
	assert_float_equal(readstr_num("test", "2-3*4", NULL), -10, 0);
	assert_float_equal(readstr_num("test", "8/2/2", NULL), 2, 0);
	assert_float_equal(readstr_num("test", "2*3+4*5-1", NULL), 25, 0);
	assert_float_equal(readstr_num("test", "1e3", NULL), 1000, 0);
	assert_float_equal(readstr_num("test", "-2*-3", NULL), 6, 0);
	assert_float_equal(readstr_num("test", "2* 3", NULL), 6, 0);//space after operator is allowed
	endptr=NULL;
	assert_float_equal(readstr_num("test", "2 *3", &endptr), 2, 0);//but not before
	assert_string_equal(endptr, " *3");
	assert_float_equal(readstr_num("test", "", NULL), 0, 0);//empty: 0
	const char *blank="   ";
	assert_float_equal(readstr_num("test", blank, &endptr), 0, 0);
	assert_ptr_equal(endptr, blank+3);//end of the blank string
	assert_float_equal(readstr_num("test", NULL, NULL), 0, 0);
	//on failure, endptr is set to the input
	const char *bad="abc";
	endptr=NULL;
	assert_true(isnan(readstr_num("test", bad, &endptr)));
	assert_ptr_equal(endptr, bad);
	const char *bad2="2*";
	endptr=NULL;
	assert_true(isnan(readstr_num("test", bad2, &endptr)));//trailing operator without number is an error
	assert_ptr_equal(endptr, bad2);
	assert_true(isnan(readstr_num("test", "1+", NULL)));
	assert_true(isnan(readstr_num("test", "1+2-", NULL)));
}
/** readstr_numarr() reads all supported types, transposes, repeats, and checks the format. */
static void readstr_numarr_types(void **state){
	(void)state;
	int nrow=0, ncol=0;
	int *ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 0, 0, M_INT, "test", "[1 2 3; 4 5 6]'"), 6);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 3);
	const int expect_t[]={1, 4, 2, 5, 3, 6};
	assert_memory_equal(ri, expect_t, sizeof(expect_t));
	free(ri); ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 0, 0, M_INT, "test", "[1.5, 2,3]"), 3);//1.5 is truncated with a warning
	assert_int_equal(ri[0], 1);
	assert_int_equal(ri[2], 3);
	free(ri); ri=NULL;
	long *rl=NULL;
	assert_int_equal(readstr_numarr((void **)&rl, &nrow, &ncol, 0, 0, M_LONG, "test", "[1 2; 3 4]'"), 4);
	const long expect_l[]={1, 3, 2, 4};
	assert_memory_equal(rl, expect_l, sizeof(expect_l));
	free(rl); rl=NULL;
	assert_int_equal(readstr_numarr((void **)&rl, &nrow, &ncol, 0, 0, M_LONG, "test", "[10000000000 2.5]"), 2);
	assert_true(rl[0]==10000000000L);
	assert_true(rl[1]==2);
	free(rl); rl=NULL;
	float *rf=NULL;
	assert_int_equal(readstr_numarr((void **)&rf, &nrow, &ncol, 0, 0, M_FLT, "test", "[0.5 1.5]'"), 2);
	assert_int_equal(nrow, 1);
	assert_int_equal(ncol, 2);
	assert_float_equal(rf[1], 1.5, 0);
	free(rf); rf=NULL;
	double *rd=NULL;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[1 2; 3 4]''"), 4);//transposed twice
	assert_int_equal(nrow, 2);
	assert_float_equal(rd[1], 2, 0);
	free(rd); rd=NULL;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[1 2 3]'*2+1"), 3);
	assert_int_equal(nrow, 1);
	assert_int_equal(ncol, 3);
	assert_float_equal(rd[2], 7, 0);
	free(rd); rd=NULL;
	//repeat with @
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[1@3 2]"), 4);
	const double expect_r[]={1, 1, 1, 2};
	assert_memory_equal(rd, expect_r, sizeof(expect_r));
	free(rd); rd=NULL;
	//automatic growth
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[0.5@25 3]"), 26);
	assert_int_equal(nrow, 26);
	assert_float_equal(rd[24], 0.5, 0);
	assert_float_equal(rd[25], 3, 0);
	free(rd); rd=NULL;
	//trailing ; is allowed
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[1 2; 3 4;]"), 4);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 2);
	free(rd); rd=NULL;
	//; before [ ends the array
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "1 2; [3]"), 2);
	assert_float_equal(rd[1], 2, 0);
	free(rd); rd=NULL;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "1 2\n3"), 2);//newline also ends it
	free(rd); rd=NULL;
	//factors before []
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "2*3/[1 2]"), 2);
	assert_float_equal(rd[0], 6, 0);
	assert_float_equal(rd[1], 3, 0);
	free(rd); rd=NULL;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", "[1 2]/2-1"), 2);
	assert_float_equal(rd[0], -0.5, 0);
	free(rd); rd=NULL;
	//empty input
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", ""), 0);
	free(rd); rd=NULL;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", NULL), 0);
	free(rd); rd=NULL;
	assert_int_equal(readstr_numarr(NULL, &nrow, &ncol, 0, 0, M_DBL, "test", "[1]"), 0);
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, 12345, "test", "[1]"), -1);//invalid type
	assert_null(rd);
}
/** readstr_numarr() with a fixed length fills, truncates, or rejects as requested. */
static void readstr_numarr_len(void **state){
	(void)state;
	int nrow=0, ncol=0;
	double *rd=(double *)malloc(5*sizeof(double));//pre-allocated
	for(int i=0; i<5; i++) rd[i]=-1;
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 5, 1, M_DBL, "test", "[3]"), 5);
	for(int i=0; i<5; i++) assert_float_equal(rd[i], 3, 0);
	assert_int_equal(nrow, 5);
	assert_int_equal(ncol, 1);
	//only the first len values are kept
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 2, 0, M_DBL, "test", "[1 2 3 4]"), 2);
	assert_float_equal(rd[0], 1, 0);
	assert_float_equal(rd[1], 2, 0);
	assert_int_equal(nrow, 2);
	//nothing supplied: kept at zero
	assert_int_equal(readstr_numarr((void **)&rd, &nrow, &ncol, 3, 0, M_DBL, "test", "[]"), 0);
	assert_int_equal(nrow, 3);
	assert_float_equal(rd[0], 0, 0);
	assert_float_equal(rd[2], 0, 0);
	free(rd); rd=NULL;
	//matrix with the right number of values
	int *ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 4, 0, M_INT, "test", "[1 2; 3 4]"), 4);
	assert_int_equal(nrow, 4);//returned as a vector when len is set
	assert_int_equal(ncol, 1);
	assert_int_equal(ri[3], 4);
	free(ri); ri=NULL;
	//not enough values
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 4, 0, M_INT, "test", "[1 2]"), -1);
	free(ri); ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 4, 1, M_INT, "test", "[1 2]"), -1);//relax=1 needs a single value
	free(ri); ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 4, 2, M_INT, "test", "[1 2]'"), -1);//cannot fill transposed
	free(ri); ri=NULL;
	assert_int_equal(readstr_numarr((void **)&ri, &nrow, &ncol, 2, 0, M_INT, "test", "[1 2 3]'"), -1);//too many to transpose
	free(ri); ri=NULL;
}
/** readstr_numarr() rejects malformed input. */
static void readstr_numarr_invalid(void **state){
	(void)state;
	const char *bad[]={
		"3 [1 2]",   //only * or / before [
		"x[1 2]",    //not a number before [
		"[1 2] abc", //garbage after ]
		"[1 2",      //unmatched [
		"1 2]",      //unmatched ]
		"[1 2; 3]",  //rows of different length
		"[1 2; 3 4; 5]",
		"[1@]",      //no count after @
		"[1 2x]",    //garbage after a number
		"1 abc",     //invalid number without []. Used to loop forever.
		"[1 abc]",   //invalid number
		"[1 2*]",    //trailing operator after a number
		"[1 2]+",    //trailing operator after ]
		"[1 2]*",
	};
	for(unsigned int i=0; i<sizeof(bad)/sizeof(bad[0]); i++){
		double *rd=NULL;
		int nrow=0, ncol=0;
		int ans=readstr_numarr((void **)&rd, &nrow, &ncol, 0, 0, M_DBL, "test", bad[i]);
		if(ans!=-1){
			fail_msg("{%s} returns %d instead of -1", bad[i], ans);
		}
		free(rd);
	}
}
/** readstr_strarr() splits by comma, honors quotes, and fills or grows the array. */
static void readstr_strarr_test(void **state){
	(void)state;
	char **res=NULL;
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "[a, b ,c]"), 3);
	assert_string_equal(res[0], "a");
	assert_string_equal(res[1], "b");
	assert_string_equal(res[2], "c");
	free_strarr(res, 3); res=NULL;
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", " ['a b', \"c,d\"] "), 2);
	assert_string_equal(res[0], "a b");
	assert_string_equal(res[1], "c,d");
	free_strarr(res, 2); res=NULL;
	//empty entries are NULL. trailing comma adds an entry.
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "a,,b,"), 4);
	assert_string_equal(res[0], "a");
	assert_null(res[1]);
	assert_string_equal(res[2], "b");
	assert_null(res[3]);
	free_strarr(res, 4); res=NULL;
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "['x',]"), 2);
	assert_string_equal(res[0], "x");
	assert_null(res[1]);
	free_strarr(res, 2); res=NULL;
	//grows beyond the initial allocation
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "1,2,3,4,5,6,7"), 7);
	assert_string_equal(res[6], "7");
	free_strarr(res, 7); res=NULL;
	//nothing
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "[]"), 0);
	assert_null(res);
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", NULL), 0);
	free(res); res=NULL;
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "  "), 0);//blank string
	assert_null(res);
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "[a,\t,b]"), 3);//blank entry
	assert_null(res[1]);
	free_strarr(res, 3); res=NULL;
	//fixed length: fill with the last entry
	assert_int_equal(readstr_strarr(&res, 3, 2, "test", "[x, y]"), 3);
	assert_string_equal(res[0], "x");
	assert_string_equal(res[1], "y");
	assert_string_equal(res[2], "y");
	for(int i=0; i<3; i++){
		free(res[i]);
	}
	assert_int_equal(readstr_strarr(&res, 3, 2, "test", "[,]"), 3);//pre-allocated. Fill NULL.
	assert_null(res[0]);
	assert_null(res[2]);
	free_strarr(res, 3); res=NULL;
	assert_int_equal(readstr_strarr(&res, 3, 1, "test", "z"), 3);
	assert_string_equal(res[2], "z");
	free_strarr(res, 3); res=NULL;
	assert_int_equal(readstr_strarr(&res, 3, 0, "test", "[x, y]"), 2);//not enough: error
	free_strarr(res, 3); res=NULL;
	//too many entries for a fixed length: error, without overflowing the array.
	assert_int_equal(readstr_strarr(&res, 2, 0, "test", "[x, y, z, w]"), 2);
	assert_string_equal(res[0], "x");
	assert_string_equal(res[1], "y");
	free_strarr(res, 2); res=NULL;
	//unmatched quote: error, keep entries before it.
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "[x, 'y]"), 1);
	assert_string_equal(res[0], "x");
	free_strarr(res, 1); res=NULL;
	assert_int_equal(readstr_strarr(&res, 3, 2, "test", "[x, \"y]"), 1);//no fill after error
	assert_string_equal(res[0], "x");
	assert_null(res[1]);
	free_strarr(res, 3); res=NULL;
	//unmatched []
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "[x, y"), 0);
	free(res); res=NULL;
	assert_int_equal(readstr_strarr(&res, 0, 0, "test", "["), 0);
	free(res); res=NULL;
}
/** trim_string() removes spaces and matching quotes without modifying the string. */
static void readstr_trim(void **state){
	(void)state;
	const char *str="  \"'abc'\"  ";
	const char *end=NULL;
	const char *start=trim_string(str, &end);
	assert_int_equal(end-start, 3);
	assert_memory_equal(start, "abc", 3);
	str="  abc  def";
	end=str+7;//limit the range
	start=trim_string(str, &end);
	assert_ptr_equal(start, str+2);
	assert_ptr_equal(end, str+5);
	end=NULL;
	assert_null(trim_string("   ", &end));
	assert_null(end);
	assert_null(trim_string(NULL, &end));
	start=trim_string(" ' ", NULL);//a single quote is kept
	assert_int_equal(start[0], '\'');
	str="'abc";
	end=NULL;
	start=trim_string(str, &end);//unmatched quote is an error and is skipped
	assert_ptr_equal(start, str+1);
	assert_ptr_equal(end, str+4);
}
/** search_keyword() finds keys preceded by a separator, ignoring case. */
static void readstr_keyword(void **state){
	(void)state;
	const char *kw=" dx=1/64 dy = 2;Dz=3,name='x y'\nscale=2*3 ";
	assert_string_equal(search_keyword(kw, "dz"), "3,name='x y'\nscale=2*3 ");
	assert_string_equal(search_keyword(kw, "DY"), "2;Dz=3,name='x y'\nscale=2*3 ");
	assert_null(search_keyword(kw, "y"));//not a key
	assert_null(search_keyword(kw, "d"));//prefix of a key
	assert_null(search_keyword(kw, "z"));
	assert_null(search_keyword(NULL, "dx"));
	assert_null(search_keyword("   ", "dx"));
	assert_float_equal(search_keyword_num(kw, "dx"), 1./64, 0);
	assert_float_equal(search_keyword_num(kw, "dy"), 2, 0);
	assert_float_equal(search_keyword_num(kw, "scale"), 6, 0);
	assert_true(isnan(search_keyword_num(kw, "dw")));
	assert_true(isnan(search_keyword_num(NULL, "dx")));
	assert_float_equal(search_keyword_num_default(kw, "dw", 5), 5, 0);
	assert_float_equal(search_keyword_num_default(kw, "dz", 5), 3, 0);
	assert_float_equal(search_keyword_num_valid(kw, "dz"), 3, 0);
	assert_true(isnan(search_keyword_num_valid(kw, "name")));//error
}

UNITTEST_MAIN(
	cmocka_unit_test(readstr_basic),
	cmocka_unit_test(readstr_array),
	cmocka_unit_test(readstr_num_edge),
	cmocka_unit_test(readstr_numarr_types),
	cmocka_unit_test(readstr_numarr_len),
	cmocka_unit_test(readstr_numarr_invalid),
	cmocka_unit_test(readstr_strarr_test),
	cmocka_unit_test(readstr_trim),
	cmocka_unit_test(readstr_keyword),
)
