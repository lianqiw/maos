/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

//To be included by read_config.c
/**
   Read in an array of double or int numbers. Supports *,/ operations
*/
int TYPEFUN1(TYPE **ret, const char *format,...)
{
    format2key;
    *ret=NULL;//initialize
    /*arrays should start with [ and end with ] or empty.*/
    TYPE data, data2;
    char *endptr,*startptr;
    int nmax;//temporary max number of elements.
    int count=0;//actual number
    nmax=10;
    long irecord=getrecord(key, 1);
    if(irecord!=-1){
	startptr=store[irecord].data;
	if(!startptr){//empty
	    return 0;
	}
	if(!(*ret=malloc(sizeof(TYPE)*nmax))){
	    error("Failed to allocate memory for ret\n");
	}

	if(startptr[0]!='['){
	    error("Invalid entry to parse for numerical array: (%s)\n", startptr);
	}
	startptr++;
	while(startptr[0]!=']' && startptr[0]!='\0'){
	    data=TYPECFUN(startptr, &endptr);
	    if(startptr==endptr){
#if STRICT==1
		info("startptr=%s, startptr[0]=%c",startptr,startptr[0]);
		error("Format is wrong:'%s'\n", (char*)store[irecord].data);
#endif
		startptr++;
	    }else{
		//no space allowed for * or /
		//while(isspace(endptr[0])) endptr++;
		while(endptr[0]=='/' || endptr[0]=='*'){
		    if(endptr[0]=='/'){
			startptr=endptr+1;
			data2=TYPECFUN(startptr, &endptr);
			if(startptr==endptr){
			    error("Error parsing '%s' for '%s' after divide /\n",startptr,TYPENAME);
			}
			data/=data2;
			//while(isspace(endptr[0])) endptr++;
		    }
		    if(endptr[0]=='*'){
			startptr=endptr+1;
			data2=TYPECFUN(startptr, &endptr);
			if(startptr==endptr){
			    error("Error parsing '%s' for '%s' after multiply *\n", startptr,TYPENAME);
			}
			data*=data2;
			//while(isspace(endptr[0])) endptr++;
		    }
		}
		(*ret)[count]=(TYPE)data;
		count++;
		if(count>=nmax){
		    nmax*=2;
		    *ret=realloc(*ret, sizeof(TYPE)*nmax);
		}
		while(endptr[0]==','||endptr[0]==';'||!isgraph((int)endptr[0])){
		    endptr++;
		}
		startptr=endptr;
	    }
	}
    }else{
	error("key '%s' not found\n", key);
	count=0;
    }
    *ret=realloc(*ret, sizeof(TYPE)*count);
    return count;
}
/**
   Read in a single double or int number. Supports *,/ operations.
*/
TYPE TYPEFUN2(const char*format,...){
    format2key;
    TYPE data, data2;
    char *endptr, *startptr;
    long irecord=getrecord(key, 1);
    if(irecord!=-1){
	startptr=store[irecord].data;
	if(!startptr){
	    error("We expect a %s, but nothing is given\n", TYPENAME);
	}
	data=TYPECFUN(startptr, &endptr);
	if(startptr==endptr){
	    error("Error parsing line '%s' for %s\n", startptr,TYPENAME);
	}
	if(strlen(endptr)>0){
	    while(endptr[0]=='/' || endptr[0]=='*'){
		if(endptr[0]=='/'){
		    startptr=endptr+1;
		    data2=TYPECFUN(startptr, &endptr);
		    if(startptr==endptr){
			error("Error parsing line '%s' for '%s' after divide /\n", 
			      startptr,TYPENAME);
		    }
		    data/=data2;
		}else if(endptr[0]=='*'){
		    startptr=endptr+1;
		    data2=TYPECFUN(startptr, &endptr);
		    if(startptr==endptr){
			error("Error parsing line '%s' for '%s' after multiply *\n", 
			      startptr,TYPENAME);
		    }
		    data*=data2;
		}
	    }
	    while(isspace((int)endptr[0])) endptr++;
	    if(strlen(endptr)>0){
		warning("Garbage found for key \"%s\":\t\"%s\"\n", key, endptr);
	    }
	}
    }else{
	error("key '%s' not found\n", key);
	data=0;
    }
    return data;
}
