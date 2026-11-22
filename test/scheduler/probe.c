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
  Connect to the scheduler at a host with scheduler_connect(), which
  authenticates TCP connections, and probe it. Used by test_auth.py.
*/
#include <unistd.h>
#include "../../sys/sys.h"
int main(int argc, char **argv){
	if(argc<2){
		printf("Usage: %s host\n", argv[0]);
		return 1;
	}
	int sock=scheduler_connect(argv[1]);
	if(sock<0){
		printf("connect failed\n");
		return 1;
	}
	int cmd[3]={CMD_PROBE, 0, 0};
	if(stwriteintarr(sock, cmd, 2)||streadintarr(sock, cmd, 3)){
		printf("probe failed\n");
		close(sock);
		return 1;
	}
	printf("probe ok %d\n", cmd[0]);
	close(sock);
	return 0;
}
