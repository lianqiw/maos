/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "common.h"
#if HAVE_LIBZSTD
#include <zstd.h>
#include <unistd.h>
#include "zstdio.h"
#include "sockio.h"
#define ZSTD_CHUNK 65536  // 64KB

struct zstd_t {
    int fd;
    int mode; // 'r' = read, 'w' = write
    union {
		struct{//for reading
        	ZSTD_DStream *dstream;//decompression
			ZSTD_inBuffer input;//for reading file {src, size, pos}
			unsigned char inbuf[ZSTD_CHUNK];//for reading file
		};
		struct{//for writting
        	ZSTD_CStream *cstream;//compression
			ZSTD_outBuffer output;//for writing to file {dst, size, pos}
			unsigned char outbuf[ZSTD_CHUNK];//for writing to file
		};
    };
};

/* Create zstd context for existing fd */
zstd_t* zstd_open(int fd, int mode) {
    if (fd < 0 || (mode != 'r' && mode != 'w')) return NULL;
    zstd_t *zf = (zstd_t*)calloc(1, sizeof(zstd_t));
    if (!zf) return NULL;
    zf->fd = fd;
    zf->mode = mode;
	
    if (mode == 'w') {
        zf->cstream = ZSTD_createCStream();
        if (!zf->cstream) { free(zf); return NULL; }
        ZSTD_initCStream(zf->cstream, 3); // default level
		zf->output.dst=zf->outbuf;
		zf->output.size=ZSTD_CHUNK;
    } else {
        zf->dstream = ZSTD_createDStream();
        if (!zf->dstream) { free(zf); return NULL; }
        ZSTD_initDStream(zf->dstream);
		zf->input.src=zf->inbuf;
    }
    return zf;
}
/*Returns 0 if success, -1 if error*/
static int handle_output(zstd_t *zf, int mode, const void *buf, size_t size) {
	if(!zf || zf->mode!='w'){
		warning("handle_output can only handle output\n");
		return -1;
	}
	ZSTD_inBuffer input = { buf, size, 0 };
	size_t ret=1;
	while((mode==0 && input.pos < input.size) || (mode!=0 && ret>0)){
		zf->output.pos=0;
		if(mode==0){
			ret = ZSTD_compressStream(zf->cstream, &zf->output, &input);
		}else if(mode==1){
			ret = ZSTD_flushStream(zf->cstream, &zf->output);
		}else{
			ret = ZSTD_endStream(zf->cstream, &zf->output);
		}
		if (ZSTD_isError(ret)) {
			warning("ZSTD returns error\n");
			return -1;//indicate error
		}
		if(zf->output.pos>0){
			if(stwrite(zf->fd, zf->output.dst, zf->output.pos)){
				warning("stwrite failed\n");
				return -1;
			}
		}
	}
	return 0;
}
/* Write buffer. Returns 0 if success, -1 if error, -2 if EOF*/
int zstd_write(zstd_t *zf, const void *buf, size_t size) {
	return handle_output(zf, 0, buf, size);
}
/* Flush buffer to file. Returns 0 if success, -1 if error,  -2 if EOF*/
int zstd_flush(zstd_t *zf) {
	return handle_output(zf, 1, 0, 0);
}
/* Read buffer.  Returns 0 if success, -1 if error,  -2 if EOF.*/
size_t zstd_read(zstd_t *zf, void *buf, size_t size) {
    if (!zf || zf->mode != 'r'){
		warning("zstd_read has wrong parameters\n");
		return 0;
	}

    ZSTD_outBuffer output = { buf, size, 0 };
    while (output.pos < output.size) {
		if(zf->input.pos >= zf->input.size){//there is no existing buffer to read
			ssize_t n = read(zf->fd, zf->inbuf, ZSTD_CHUNK);
			if (n == 0){
				//warning("zstd_read: end of file for size=%zu\n", size);
				break; // EOF
			}
			if (n < 0) {
				if (errno == EINTR) continue;
				warning("zstd_read: reading error for size=%zu\n", size);
				return -1;
			}
			zf->input.pos=0;
			zf->input.size=n;
		}
        
        size_t ret = ZSTD_decompressStream(zf->dstream, &output, &zf->input);
        if (ZSTD_isError(ret)) {
			warning("ZSTD_decompressStream returns error\n");
			return -1;
		}
    }
    return convert_ans(output.pos, size);
}



/* Close wrapper */
int zstd_close(zstd_t *zf, int close_fd) {
    if (!zf) return -1;
    int err = 0;

    if (zf->mode == 'w') {
		handle_output(zf, 2, 0, 0);
        ZSTD_freeCStream(zf->cstream);
    } else {
        ZSTD_freeDStream(zf->dstream);
    }

    if (close_fd) close(zf->fd);
    free(zf);
    return err;
}

#endif
