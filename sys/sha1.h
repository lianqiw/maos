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
#ifndef AOS_SYS_SHA1_H
#define AOS_SYS_SHA1_H
typedef struct SHA1_CTX SHA1_CTX;
void sha1_init(SHA1_CTX *ctx);
void sha1_update(SHA1_CTX *ctx, const unsigned char *data, size_t len);
void sha1_final(SHA1_CTX *ctx, unsigned char *out);
void base64_encode(const unsigned char *in, int in_len, char *out);
void base64_sha1(const char* in,  char *out);
#endif
