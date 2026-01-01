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
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
/* Minimal SHA-1 Implementation */
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "sha1.h"
struct SHA1_CTX{
    uint32_t h[5];
    unsigned char block[64];
    uint64_t bitlen;
    size_t curlen;
};

void sha1_init(SHA1_CTX *ctx) {
    ctx->h[0] = 0x67452301;
    ctx->h[1] = 0xEFCDAB89;
    ctx->h[2] = 0x98BADCFE;
    ctx->h[3] = 0x10325476;
    ctx->h[4] = 0xC3D2E1F0;
    ctx->curlen = 0;
    ctx->bitlen = 0;
}

static void sha1_transform(SHA1_CTX *ctx, const unsigned char *data) {
    uint32_t w[80];
    uint32_t a, b, c, d, e, t;

    for (int i = 0; i < 16; i++)
        w[i] = (data[4*i]<<24)|(data[4*i+1]<<16)|(data[4*i+2]<<8)|data[4*i+3];
    for (int i = 16; i < 80; i++) {
        uint32_t v = w[i-3]^w[i-8]^w[i-14]^w[i-16];
        w[i] = (v << 1) | (v >> 31);
    }

    a = ctx->h[0]; b = ctx->h[1]; c = ctx->h[2]; d = ctx->h[3]; e = ctx->h[4];

    for (int i = 0; i < 80; i++) {
        uint32_t f, k;
        if (i < 20) { f = (b & c) | ((~b) & d); k = 0x5A827999; }
        else if (i < 40) { f = b ^ c ^ d; k = 0x6ED9EBA1; }
        else if (i < 60) { f = (b & c) | (b & d) | (c & d); k = 0x8F1BBCDC; }
        else { f = b ^ c ^ d; k = 0xCA62C1D6; }

        t = ((a << 5) | (a >> 27)) + f + e + k + w[i];
        e = d; d = c; c = (b << 30) | (b >> 2); b = a; a = t;
    }

    ctx->h[0] += a; ctx->h[1] += b; ctx->h[2] += c;
    ctx->h[3] += d; ctx->h[4] += e;
}

void sha1_update(SHA1_CTX *ctx, const unsigned char *data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        ctx->block[ctx->curlen++] = data[i];
        if (ctx->curlen == 64) {
            sha1_transform(ctx, ctx->block);
            ctx->bitlen += 512;
            ctx->curlen = 0;
        }
    }
}

void sha1_final(SHA1_CTX *ctx, unsigned char *out) {
    ctx->bitlen += ctx->curlen * 8;
    ctx->block[ctx->curlen++] = 0x80;
    if (ctx->curlen > 56) {
        while (ctx->curlen < 64) ctx->block[ctx->curlen++] = 0;
        sha1_transform(ctx, ctx->block);
        ctx->curlen = 0;
    }
    while (ctx->curlen < 56) ctx->block[ctx->curlen++] = 0;

    for (int i = 7; i >= 0; i--)
        ctx->block[ctx->curlen++] = (ctx->bitlen >> (i * 8)) & 0xFF;

    sha1_transform(ctx, ctx->block);
    for (int i = 0; i < 5; i++) {
        out[i*4] = (ctx->h[i] >> 24) & 0xFF;
        out[i*4+1] = (ctx->h[i] >> 16) & 0xFF;
        out[i*4+2] = (ctx->h[i] >> 8) & 0xFF;
        out[i*4+3] = ctx->h[i] & 0xFF;
    }
}

/**  Base64 Encode  */
void base64_encode(const unsigned char *in, int in_len, char *out) {
	static const char b64_table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int i, j;
    for (i = 0, j = 0; i < in_len;) {
        uint32_t a = i < in_len ? in[i] : 0; i++;
        uint32_t b = i < in_len ? in[i] : 0; i++;
        uint32_t c = i < in_len ? in[i] : 0; i++;
        uint32_t triple = (a << 16) | (b << 8) | c;

        out[j++] = b64_table[(triple >> 18) & 0x3F];
        out[j++] = b64_table[(triple >> 12) & 0x3F];
        out[j++] = (i > in_len + 1) ? '=' : b64_table[(triple >> 6) & 0x3F];
        out[j++] = (i > in_len) ? '=' : b64_table[triple & 0x3F];
    }
    out[j] = 0;
}
void base64_sha1(const char* in,  char *out){
	unsigned char sha[20];
	SHA1_CTX sha_ctx;
    sha1_init(&sha_ctx);
    sha1_update(&sha_ctx, (unsigned char*)in, strlen(in));
    sha1_final(&sha_ctx, sha);
    base64_encode(sha, sizeof(sha), out);
}
