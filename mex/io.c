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
#include <stdarg.h>
#include <ctype.h>
#include"io.h"
static void write_timestamp(file_t *fp);

static const char *myasctime(void){
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
    if(!check_suffix(fn2,".bin") && !check_suffix(fn2, ".bin.gz")
       &&!check_suffix(fn2,".fits") && !check_suffix(fn2, ".fits.gz")){
	strncat(fn2, ".bin", 4);
    }
    if(mod[0]=='r' || mod[0]=='a'){
	if(!exist(fn2)){/*If does not exist.*/
	    if(!check_suffix(fn2, ".gz")){
		strncat(fn2, ".gz", 3);
	    }else if (check_suffix(fn, ".gz")){
		/*ended with bin.gz, change to .bin*/
		fn2[strlen(fn2)-3]='\0';
	    }
	    if(!exist(fn2)){
		return NULL;
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
    if((check_suffix(fn, ".bin") || check_suffix(fn, ".fits")) && mod[0]=='w'){
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
    if(check_suffix(fn, ".fits") || check_suffix(fn, ".fits.gz")){
	fp->isfits=1;
    }
    if(mod[0]=='w' && !fp->isfits){
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
static inline void zfwrite_do(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    if(fp->isgzip){
	if(gzwrite((voidp)fp->p, ptr, size*nmemb)!=size*nmemb){
	    perror("gzwrite");
	    error("write failed\n");
	}
    }else{
	if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    perror("fwrite");
	    error("writefailed\n");
	}
    }
}
/**
   Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
   fwrite. Follows the interface of fwrite.
 */
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits && size>1){
	/* write a block of 2880 bytes each time, with big-endianness.*/
	const int bs=2880;
	char junk[bs];
	int length=size*nmemb;
	int nb=(length+bs-1)/bs;
	char *in=(char*)ptr;
	for(int ib=0; ib<nb; ib++){
	    int nd=length<bs?length:bs;
	    switch(size){
	    case 2:
		for(int i=0; i<nd; i+=2){
		    junk[i]=in[i+1];
		    junk[i+1]=in[i];
		}
		break;
	    case 4:
		for(int i=0; i<nd; i+=4){
		    junk[i]  =in[i+3];
		    junk[i+1]=in[i+2];
		    junk[i+2]=in[i+1];
		    junk[i+3]=in[i  ];
		}
		break;
	    case 8:
	    case 16:
		for(int i=0; i<nd; i+=8){
		    junk[i  ]=in[i+7];
		    junk[i+1]=in[i+6];
		    junk[i+2]=in[i+5];
		    junk[i+3]=in[i+4];
		    junk[i+4]=in[i+3];
		    junk[i+5]=in[i+2];
		    junk[i+6]=in[i+1];
		    junk[i+7]=in[i  ];
		}
		break;
	    default:
		error("Invalid\n");
	    }
	    /* use bs instead of nd to test tailing blanks*/
	    in+=bs; length-=bs;
	    if(length<0){
		memset(junk+nd, 0, (bs-nd)*sizeof(char));
	    }
	    zfwrite_do(junk, sizeof(char), bs, fp);
	}
    }else{
	zfwrite_do(ptr, size, nmemb, fp);
    }
}
void zfwrite_complex(const double* pr, const double *pi,const size_t nmemb, file_t *fp){
    dcomplex *tmp=malloc(sizeof(dcomplex)*nmemb);
    long i;
    for(i=0; i<nmemb; i++){
	tmp[i]=pr[i]+I*pi[i];
    }
    zfwrite(tmp, sizeof(dcomplex), nmemb, fp);
    free(tmp);
}

static inline int zfread_do(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(fp->isgzip){
	return gzread((voidp)fp->p, ptr, size*nmemb)>0?0:-1;
    }else{
	return fread(ptr, size, nmemb, (FILE*)fp->p)==nmemb?0:-1;
    }
}
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread. It does byte ordering from big endian
   to small endian in case we are reading fits file.  */
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits && size>1){/*need to do byte swapping.*/
	const long bs=2880;
	char junk[bs];
	long length=size*nmemb;
	long nb=(length+bs-1)/bs;
	char *out=(char*)ptr;
	for(int ib=0; ib<nb; ib++){
	    if(zfread_do(junk, sizeof(char), bs, fp)) return -1;
	    int nd=length<bs?length:bs;
	    switch(size){
	    case 2:
		for(int i=0; i<nd; i+=2){
		    out[i]  =junk[i+1];
		    out[i+1]=junk[i];
		}
		break;
	    case 4:
		for(int i=0; i<nd; i+=4){
		    out[i]  =junk[i+3];
		    out[i+1]=junk[i+2];
		    out[i+2]=junk[i+1];
		    out[i+3]=junk[i  ];
		}
		break;
	    case 8:
	    case 16:
		for(int i=0; i<nd; i+=8){
		    out[i  ]=junk[i+7];
		    out[i+1]=junk[i+6];
		    out[i+2]=junk[i+5];
		    out[i+3]=junk[i+4];
		    out[i+4]=junk[i+3];
		    out[i+5]=junk[i+2];
		    out[i+6]=junk[i+1];
		    out[i+7]=junk[i  ];
		}
		break;
	    default:
		error("Invalid\n");
	    }
	    out+=bs;
	    length-=bs;
	}
	return 0;
    }else{
	return zfread_do(ptr, size, nmemb, fp);
    }
}
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread.
 */
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(zfread2(ptr, size, nmemb, fp)){
	error("Error happened while reading\n");
    }
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

int zfeof(file_t *fp){
    return zfseek(fp, 1, SEEK_SET)<0?-1:0;
}
/**
   Write the magic into file. Also write a dummy header to make data alignment to 8 bytes.
*/
static void write_bin_magic(uint32_t magic, file_t *fp){
    if(fp->isfits) error("fits file is not supported\n");
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
static void write_bin_header(const char *header, file_t *fp){
    if(!header) return;
    if(fp->isfits) error("fits file is not supported\n");
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

static void write_timestamp(file_t *fp){
    if(fp->isfits) error("Not supported\n");
    char header[128];
    snprintf(header,128, "Created by write on %s\n",
	     myasctime());
    write_bin_header(header, fp);
}

/**
   Obtain the current magic number. If it is a header, read it out if output of
header is not NULL.  The header will be appended to the output header.*/
static uint32_t read_bin_magic(file_t *fp, char **header){
    uint32_t magic,magic2;
    uint64_t nlen, nlen2;
    if(fp->isfits) error("fits file is not supported\n");
    while(1){
	/*read the magic number.*/
	if(zfread2(&magic, sizeof(uint32_t), 1, fp)) return 0;
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
	      error("Header verification failed: %u<>%u, %lu<>%lu\n", magic, magic2, (unsigned long)nlen, (unsigned long)nlen2);
	    }
	}else{ /*otherwise return the magic number*/
	    return magic;
	}
    }/*while*/
}
/**
   Write fits header. extra is extra header that will be put in fits comment
*/
static void
write_fits_header(file_t *fp, const char *str, uint32_t magic, int count, ...){
    uint64_t naxis[count];
    va_list ap;
    va_start (ap, count);              /*Initialize the argument list. */
    int empty=0;
    for (int i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	if((*addr)==0) empty=1;
	naxis[i]=*addr;
    }
    va_end(ap);

    if(empty) count=0;

    int bitpix;
    switch(magic){
    case M_FLT:
	bitpix=-32;
	break;
    case M_DBL:
	bitpix=-64;
	break;
    default:
	error("Data type is not yet supported.\n");
    }
    const int nh=36;
    char header[nh][80];
    memset(header, ' ', sizeof(char)*36*80);
    int hc=0;
    if(fp->isfits==1){
	snprintf(header[hc], 80, "%-8s= %20s", "SIMPLE", "T");    header[hc][30]=' '; hc++;
	fp->isfits++;
    }else{
	snprintf(header[hc], 80, "%-8s= %s", "XTENSION", "'IMAGE   '");    header[hc][20]=' '; hc++;
    }
    snprintf(header[hc], 80, "%-8s= %20d", "BITPIX", bitpix); header[hc][30]=' '; hc++;
    snprintf(header[hc], 80, "%-8s= %20d", "NAXIS", count);   header[hc][30]=' '; hc++;
#define FLUSH_OUT /*write the block and reset */	\
	if(hc==nh){					\
	    zfwrite(header, sizeof(char), 36*80, fp);	\
	    memset(header, ' ', sizeof(char)*36*80);	\
	    hc=0;					\
	}
    for (int i = 0; i < count; i++){
	FLUSH_OUT;
	snprintf(header[hc], 80, "%-5s%-3d= %20lu", "NAXIS", i+1, (unsigned long)(naxis[i])); header[hc][30]=' '; hc++;
    }
    if(str){
	const char *str2=str+strlen(str);
	while(str<str2){
	    char *nl=strchr(str, '\n');
	    int length;
	    if(nl){
		length=nl-str+1;
		if(length<=70) nl[0]=';';
	    }else{
		length=strlen(str);
	    }
	    if(length>70) length=70;
	    FLUSH_OUT;
	    strncpy(header[hc], "COMMENT   ", 10);
	    strncpy(header[hc]+10, str, length);
	    hc++;
	    str+=length;
	}
    }
    FLUSH_OUT;
    snprintf(header[hc], 80, "%-8s", "END"); header[hc][8]=' '; hc++;
    zfwrite(header, sizeof(char), 36*80, fp);
#undef FLUSH_OUT
}

/**
   Read fits header
 */
int read_fits_header(file_t *fp, char **str, uint32_t *magic, uint64_t *nx, uint64_t *ny){
    char line[81];
    int end=0;
    int page=0;
    int bitpix=0;
    int naxis=0;
    while(!end){
	int start=0;
	if(page==0){
	    if(zfread2(line, 1, 80, fp)) return -1; line[80]='\0';
	    if(strncmp(line, "SIMPLE", 6) && strncmp(line, "XTENSION= 'IMAGE", 16)){
		error("Fits header is not recognized\n");
	    }
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &bitpix)!=1) error("Unable to determine bitpix\n");
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &naxis)!=1) error("Unable to determine naxis\n");
	    if(naxis>2) error("Data type not supported\n");
	    if(naxis>0){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", nx)!=1) error("Unable to determine nx\n");
	    }else{
		*nx=0;
	    }
	    if(naxis>1){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", ny)!=1) error("Unable to determine ny\n");
	    }else{
		*ny=0;
	    }
	    start=3+naxis;
	}
	for(int i=start; i<36; i++){
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(!strncmp(line, "COMMENT", 7)){
		for(int j=79; j>9; j--){
		    if(isspace((int)line[j])){
			line[j]='\0';
		    }else{
			break;
		    }
		}
		int length=strlen(line+10);
		if(*str){
		    *str=realloc(*str, strlen(*str)+length+1);
		}else{
		    *str=malloc(length+1); (*str)[0]='\0';
		}
		strcat(*str, line+10);
	    }else if(!strncmp(line, "END",3)){
		end=1;
	    }
	}
    }
    switch(bitpix){
    case -32:
	*magic=M_FLT;
	break;
    case -64:
	*magic=M_DBL;
	break;
    default:
	error("Invalid\n");
    }
    return 0;
}

/**
   Write multiple long numbers into the file. To write three numbers, a, b, c,
call with zfwritelarr(fp, 3, &a, &b, &c); */
void zfwritelarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              /*Initialize the argument list. */
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	zfwrite(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       /* Clean up.  */
}
/**
   Read multiple long numbers from the file. To read three numbers, a, b, c,
   call with zfreadlarr(fp, 3, &a, &b, &c);
 */
void zfreadlarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              /*Initialize the argument list. */
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	zfread(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       /* Clean up.  */
}
/**
  A unified header writing routine for .bin and .fits files. It write the array
information and string header if any.  */
void write_header(const header_t *header, file_t *fp){
    if(fp->isfits){
	if(!iscell(header->magic)){
	    write_fits_header(fp, header->str, header->magic, 2, &header->nx, &header->ny);
	}
    }else{
	if(header){
	    write_bin_header(header->str, fp);
	}
	write_bin_magic(header->magic, fp);
	zfwritelarr(fp, 2, &header->nx, &header->ny);
    }
}
/**
   A unified header reading routine for .bin and .fits files. It read the array
information and string header if any.  Return error signal.*/
int read_header2(header_t *header, file_t *fp){
    int ans;
    header->str=NULL;
    if(fp->isfits){
	ans=read_fits_header(fp, &header->str, &header->magic, &header->nx, &header->ny);
    }else{
	header->magic=read_bin_magic(fp, &header->str);
	if(header->magic==0){
	    ans=-1;
	}else{
	    ans=0;
	    zfreadlarr(fp, 2, &header->nx, &header->ny);
	}
    }
    return ans;
}
/**
   calls read_header2 and abort if error happens.*/
void read_header(header_t *header, file_t *fp){
    if(read_header2(header, fp)){
	error("read_header failed\n");
    }
}
