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
  Test the job list and CPU/GPU accounting of the scheduler
  (tools/scheduler_job.c) without any socket. monitor_send() is replaced by a
  stub that records the updates sent to monitors. No job is launched.
*/
#include <sys/wait.h>
#include "../tools/scheduler_job.c"
#include "unittest.h"

static int nsend=0;//number of monitor_send() calls
static int last_info=0;//status of the last job sent
static int last_pid=0;//pid of the last job sent
void monitor_send(RUN_T* irun, const char* path){
	(void)path;
	nsend++;
	last_info=irun->status.info;
	last_pid=irun->pid;
}
/** Return the pid of a process that has exited. */
static int dead_pid(void){
	int pid=fork();
	if(pid==0) _exit(0);
	waitpid(pid, NULL, 0);
	return pid;
}

/** Queue a job, kill it while queued, restart it, and remove it. */
static void job_queue(void **state){
	(void)state;
	assert_null(running_head());
	assert_null(runned_head());
	nsend=0;
	queue_new_job("/bin/true", "maos\na=1");
	RUN_T *irun=running_head();
	assert_non_null(irun);
	assert_true(irun->pid<0);
	assert_int_equal(irun->status.info, S_QUEUED);
	assert_string_equal(irun->exe, "/bin/true");
	assert_string_equal(irun->path0, "maos\na=1");
	assert_string_equal(irun->path, "maos a=1");
	assert_int_equal(irun->sock, -1);
	assert_int_equal(nsend, 2);//path and status
	assert_ptr_equal(running_get_by_pid(irun->pid), irun);
	assert_null(running_get_by_pid(0));//0 matches a launched job only
	int pid=irun->pid;

	running_kill(pid);//queued job is removed as killed
	assert_null(running_head());
	assert_ptr_equal(runned_head(), irun);
	assert_int_equal(irun->status.info, S_KILLED);
	assert_int_equal(irun->status.done, 1);
	assert_int_equal(last_info, S_KILLED);
	assert_int_equal(nrun_get(0), 0);

	runned_restart(pid);//back to the queue with a new pid
	assert_null(runned_head());
	assert_ptr_equal(running_head(), irun);
	assert_int_equal(irun->status.info, S_QUEUED);
	assert_true(irun->pid<pid);
	assert_int_equal(irun->pidnew, irun->pid);
	assert_int_equal(last_pid, pid);//monitor is told with the old pid
	pid=irun->pid;

	running_remove(pid, S_FINISH);
	assert_ptr_equal(runned_head(), irun);
	runned_remove(pid);
	assert_null(runned_head());
	runned_remove(pid);//not found: only a warning
	runned_restart(pid);
}

/** A launched job is counted once for its CPUs and GPUs and released when it ends. */
static void job_nrun(void **state){
	(void)state;
	int pid=getpid();
	RUN_T *irun=running_add(pid, 5);
	assert_non_null(irun);
	assert_non_null(irun->exe);
	assert_ptr_equal(running_add(pid, 6), irun);//only the socket is updated
	assert_int_equal(irun->sock, 6);
	assert_ptr_equal(running_get_by_sock(6), irun);
	assert_null(running_get_by_sock(5));
	assert_ptr_equal(running_get_by_pid(0), irun);
	irun->nthread=3;
	irun->ngpu=1;
	irun->status.info=S_START;
	nrun_add(irun);
	nrun_add(irun);
	assert_int_equal(nrun_get(0), 3);
	assert_int_equal(nrun_get(1), 1);
	nrun_sub(irun);
	nrun_sub(irun);
	assert_int_equal(nrun_get(0), 0);
	assert_int_equal(nrun_get(1), 0);
	nrun_add(irun);
	running_remove(pid, S_FINISH);//releases the resource
	assert_int_equal(irun->counted, 0);
	assert_int_equal(nrun_get(0), 0);
	assert_int_equal(nrun_get(1), 0);
	assert_int_equal(irun->status.info, S_FINISH);
	running_remove(pid, S_FINISH);//not found: no effect
	runned_remove(-1);
	assert_null(runned_head());
	assert_null(running_add(dead_pid(), -1));//already exited
	assert_null(running_head());
}

/** A job whose process is gone is marked crashed after the grace period. */
static void job_exists(void **state){
	(void)state;
	RUN_T *irun=running_add(getpid(), -1);
	assert_non_null(irun);
	irun->status.info=S_START;
	irun->nthread=2;
	nrun_add(irun);
	irun->last_time=myclocki();
	assert_int_equal(check_job_exists(irun), 1);
	irun->pid=irun->pidnew=dead_pid();
	assert_int_equal(check_job_exists(irun), 1);//within grace period
	assert_ptr_equal(running_head(), irun);
	irun->last_time=myclocki()-61;
	assert_int_equal(check_job_exists(irun), 0);
	assert_null(running_head());
	assert_ptr_equal(runned_head(), irun);
	assert_int_equal(irun->status.info, S_CRASH);
	assert_int_equal(nrun_get(0), 0);
	assert_int_equal(check_job_exists(irun), 0);//finished job is not checked
	runned_remove(-1);
	check_jobs();
	assert_int_equal(nrun_get(0), 0);
}

/** Newlines in the job path are replaced by spaces. */
static void job_remove_endl(void **state){
	(void)state;
	char *path=remove_endl("a\nb\n\nc");
	assert_string_equal(path, "a b  c");
	free(path);
	path=remove_endl("");
	assert_string_equal(path, "");
	free(path);
}

UNITTEST_MAIN(
	cmocka_unit_test(job_queue),
	cmocka_unit_test(job_nrun),
	cmocka_unit_test(job_exists),
	cmocka_unit_test(job_remove_endl),
)
