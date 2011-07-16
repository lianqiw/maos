#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include"io.h"
const char *myasctime(void){
    static char st[64];
    time_t a;
    time(&a);
    ctime_r(&a, st);
    st[strlen(st)-1]='\0';
    return st;
}
static int exist(const char *fn){
    /**
       Test whether a file exists.
    */
    struct stat buf;
    return !stat(fn, &buf);
}
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.
 */
static int mystrcmp(const char *a, const char *b){
    if(!a || !b) return 1;
    int la=strlen(a);
    int lb=strlen(b);
    if(la==0 || lb==0 || la<lb){
	return 1;
    }else{
	return strncmp(a,b,lb);
    }
}
/**
   Check the suffix of a file.
*/
static int check_suffix(const char *fn, const char *suffix){
    if(!fn || !suffix) return 0;
    int lfn=strlen(fn);
    int lsu=strlen(suffix);
    if(lfn < lsu) return 0;
    if(mystrcmp(fn+lfn-lsu,suffix)){
	return 0;
    }else{
	return 1;
    }
}
/**
   Test whether fn is a symbolic link
*/
static int islink(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISLNK(buf.st_mode);
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
    /*If there is no recognized suffix, add .bin in the end.*/
    if(!check_suffix(fn2,".bin") && !check_suffix(fn2, ".bin.gz")){
	strncat(fn2, ".bin", 4);
    }
    if(mod[0]=='r' || mod[0]=='a'){
	if(!exist(fn2)){/*If does not exist.*/
	    if(check_suffix(fn2, ".bin")){
		/*ended with .bin, change to .bin.gz*/
		strncat(fn2, ".gz", 3);
		if(!exist(fn2)){
		    return NULL;
		}
	    }else{
		/*ended with bin.gz, change to .bin*/
		fn2[strlen(fn2)-3]='\0';
		if(!exist(fn2)){
		    return NULL;
		}
	    }
	}
    }else if (mod[0]=='w'){
	if(islink(fn2)){
	    /*remove old file to avoid write over a symbolic link.*/
	    if(remove(fn2)){
		error("Failed to remove %s\n", fn2);
	    }
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
}
/*stripped down version of io.c*/
file_t* zfopen(const char *fn_in, char *mod){
    char *fn=procfn(fn_in, mod, 1);
    if(!fn){
	error("%s does not exist\n", fn_in);
    }
    file_t* fp=calloc(1, sizeof(file_t));
    if(check_suffix(fn, ".bin") && mod[0]=='w'){
	fp->isgzip=0;
	if(!(fp->p=fopen(fn,mod))){
	    error("Error fopen for %s\n",fn);
	}
    }else{ 
	fp->isgzip=1;
	if(!(fp->p=gzopen(fn,mod))){
	    error("Error gzopen for %s\n",fn);
	}
    }
    if(mod[0]=='w'){
	write_timestamp(fp);
    }
    free(fn);
    return fp;
}
void zfclose(file_t *fp){
    if(fp->isgzip){
	gzclose((voidp)fp->p);
    }else{
	if(fclose((FILE*)fp->p)){
	    perror("fclose\n");
	}
    }
    free(fp);
}
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isgzip){
	gzwrite((voidp)fp->p, ptr, size*nmemb);
    }else{
	if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    error("Write failed\n");
	}
    }
}
void zfwrite_complex(const double* pr, const double *pi,const size_t nmemb, file_t *fp){
    dcomplex *tmp=malloc(sizeof(dcomplex)*nmemb);
    long i;
    for(i=0; i<nmemb; i++){
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
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->eof) return;
    if(fp->isgzip){
	int status;
	if((status=gzread((voidp)fp->p, ptr, size*nmemb))<1){
	    if(status==-1){
		warning("Error happened in reading");
	    }else{
		warning("End of File encoutnered");
	    }
	    fp->eof=1;
	}
    }else{
	size_t nmemb2;
	if((nmemb2=fread(ptr, size, nmemb, (FILE*)fp->p))!=nmemb){
	    fp->eof=1;
	    if(feof((FILE*)fp->p)){
		warning("End of File encountered!. Want %lu, get %lu\n", nmemb, nmemb2);
	    }else{
		warning("Read failed. Unknown error.\n");
	    }
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
/**
   Move the current position pointer, like fseek
*/
int zfseek(file_t *fp, long offset, int whence){
    if(fp->isgzip){
	int res=gzseek((voidp)fp->p,offset,whence);
	if(res<0)
	    return 1;
	else
	    return 0;
    }else{
	return fseek((FILE*)fp->p,offset,whence);
    }
}
/**
   Write the magic into file. Also write a dummy header to make data alignment to 8 bytes.
*/
void write_magic(uint32_t magic, file_t *fp){
    uint32_t magic2=M_SKIP;
    zfwrite(&magic2, sizeof(uint32_t), 1, fp);
    zfwrite(&magic,  sizeof(uint32_t), 1, fp);
}
/**
   Append to the file the header to the end of the file(or rather, the
   tailer). First write magic number, then the length of the header, then the
   header, then the length of the header again, then the magic number again. The
   two set of strlen and header are used to identify the header from the end of
   the file and also to verify that what we are reading are indeed header. The
   header may be written multiple times. They will be concatenated when
   read. The header should contain key=value entries just like the configuration
   files. The entries should be separated by new line charactor. */
void write_header(const char *header, file_t *fp){
    if(!header) return;
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(header)+1;
    /*make header 8 byte alignment.*/
    char *header2=strdup(header);
    if(nlen % 8 != 0){
	nlen=(nlen/8+1)*8;
	header2=realloc(header2, nlen);
    }
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(header2, 1, nlen, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    free(header2);
}
void write_timestamp(file_t *fp){
    char header[128];
    snprintf(header,128, "Created by write on %s\n",
	     myasctime());
    write_header(header, fp);
}
/**
   Obtain the current magic number. If it is a header, read it out if output of
header is not NULL.  The header will be appended to the output header.*/
uint32_t read_magic(file_t *fp, char **header){
    uint32_t magic,magic2;
    uint64_t nlen, nlen2;
    while(1){
	/*read the magic number.*/
	zfread(&magic, sizeof(uint32_t), 1, fp);
	/*If it is header, read or skip it.*/
	if(magic==M_SKIP){
	    continue;
	}else if(magic==M_HEADER){
	    zfread(&nlen, sizeof(uint64_t), 1, fp);
	    if(nlen>0){
		if(header){
		    char header2[nlen];
		    zfread(header2, 1, nlen, fp);
		    header2[nlen-1]='\0'; /*make sure it is NULL terminated.*/
		    if(*header){
			*header=realloc(*header, ((*header)?strlen(*header):0)+strlen(header2)+1);
			strncat(*header, header2, nlen);
		    }else{
			*header=strdup(header2);
		    }
		}else{
		    zfseek(fp, nlen, SEEK_CUR);
		}
	    }
	    zfread(&nlen2, sizeof(uint64_t),1,fp);
	    zfread(&magic2, sizeof(uint32_t),1,fp);
	    if(magic!=magic2 || nlen!=nlen2){
		error("Header verification failed: %u<>%u, %lu<>%lu\n", magic, magic2, nlen, nlen2);
	    }
	}else{ /*otherwise return the magic number*/
	    return magic;
	}
    }/*while*/
}
