#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>

#include"io.h"

static int exist(const char *fn){
    /**
       Test whether a file exists.
    */
    struct stat buf;
    return !stat(fn, &buf);
}

static char* procfn(const char *fn, const char *mod,const int gzip){
    char *fn2;
    if(fn[0]=='~'){
	char *HOME=getenv("HOME");
	fn2=malloc(strlen(HOME)+strlen(fn)+16);
	strcpy(fn2,HOME);
	strcat(fn2,fn+1);
    }else{
	fn2=malloc(strlen(fn)+16);
	strcpy(fn2,fn);
    }
    if(mod[0]=='r'){
	if(!exist(fn2)){
	    if(strlen(fn2)>=7&&!strncmp(fn2+strlen(fn2)-7,".bin.gz",7)){
		//ended with bin.gz
		fn2[strlen(fn2)-3]='\0';
		if(!exist(fn2)){
		    error("Neither %s nor %s exist\n", fn,fn2);
		}
	    }else if(strlen(fn2)>=4&&!strncmp(fn2+strlen(fn2)-4,".bin",4)){
		//ended with .bin
		strncat(fn2, ".gz", 3);
		if(!exist(fn2)){
		    error("Neither %s nor %s exist\n", fn,fn2);
		}
	    }else{//no recognized suffix
		strncat(fn2, ".bin", 4);
		if(!exist(fn2)){
		    strncat(fn2,".gz",3);
		    if(!exist(fn2)){
			error("Neither %s, %s.bin, nor %s.bin.gz exist\n",
			      fn,fn,fn);
		    }
		}
	    }
	}
    }else if (mod[0]=='w'){//for write, no suffix. we append .bin.gz
	if(!((strlen(fn2)>=7&&!strncmp(fn2+strlen(fn2)-7,".bin.gz",7))
	     || (strlen(fn2)>=4&&!strncmp(fn2+strlen(fn2)-4,".bin",4)))){
	    if(gzip)
		strncat(fn2,".bin.gz",7);
	    else
		strncat(fn2,".bin",4);
	}
	if(exist(fn2)){//remove old file to avoid write over a symbolic link.
	    remove(fn2);
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
}
/*stripped down version of io.c*/
file_t* openfile(const char *fn_in, char *mod){
    char *fn=procfn(fn_in, mod, 1);
    file_t* fp=calloc(1, sizeof(file_t));
    if(!strcmp(fn+strlen(fn)-3,".gz")){
	fp->isgzip=1;
	if(!(fp->p=gzopen(fn,mod))){
	    error("Error gzopen for %s\n",fn);
	}
    }else{ 
	fp->isgzip=0;
	if(!(fp->p=fopen(fn,mod))){
	    error("Error fopen for %s\n",fn);
	}
    }
    free(fn);
    return fp;
}
void closefile(file_t *fp){
    if(fp->isgzip){
	gzclose((voidp)fp->p);
    }else{
	fclose((FILE*)fp->p);
    }
    free(fp);
}
void writefile(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isgzip){
	gzwrite((voidp)fp->p, ptr, size*nmemb);
    }else{
	if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    error("Write failed\n");
	}
    }
}
void writefile_complex(const double* pr, const double *pi,const size_t nmemb, file_t *fp){
    dcomplex *tmp=malloc(sizeof(dcomplex)*nmemb);
    for(long i=0; i<nmemb; i++){
	tmp[i]=pr[i]+I*pi[i];
    }
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isgzip){
	gzwrite((voidp)fp->p, tmp, sizeof(dcomplex)*nmemb);
    }else{
	if(fwrite(tmp, sizeof(dcomplex), nmemb, (FILE*)fp->p)!=nmemb){
	    error("Write failed\n");
	}
    }
    free(tmp);
}
void readfile(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->eof) return;
    if(fp->isgzip){
	int status;
	if((status=gzread((voidp)fp->p, ptr, size*nmemb))<1){
	    if(status==-1){
		warning("Error happened in reading");
	    }
	    else{
		warning("End of File encoutnered");
		fp->eof=1;
	    }
	}
    }else{
	if(fread(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    error("Read failed\n");
	}
	if(feof((FILE*)fp->p)){
	    warning("End of File encountered!\n");
	    fp->eof=1;
	}
    }
}

int test_eof(file_t *fp){
    int ans;
    long cpos, fpos;
    if(fp->isgzip){
	cpos=gztell((voidp)fp->p);
	gzseek((voidp)fp->p,0,SEEK_END);
	fpos=gztell((voidp)fp->p);
    }else{
	cpos=ftell((FILE*)fp->p);
	fseek((FILE*)fp->p, 0, SEEK_END);    
	fpos=ftell((FILE*)fp->p);
    }
    if(cpos!=fpos){
	fprintf(stderr,"There are unread bytes near the end of the file\n"
		"Current position is %ld, File END is %ld\n",cpos, fpos);
	ans=1;
    }else{
	ans=0;
    }
    return ans;
}
