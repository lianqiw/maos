\page page94_scheduler Scheduler Protocol

This page describes how \c maos, \c monitor, \c drawdaemon, the browser, and
\c scripts/maos_client.py talk to the \c scheduler. See \ref compile-browser
and \ref compile-monitor for how to use them. The command values are defined
in \c sys/scheduler_client.h. Existing values must not be changed, because
\c maos and \c monitor of different versions connect to the same scheduler.

\section sched-files Source Files

- \c tools/scheduler.c &ndash; \c main(), the command handlers
  (\c respond_* listed in \c respond_list), and the monitor list.
- \c tools/scheduler_job.c &ndash; the list of queued, running, and
  finished jobs and the CPU/GPU accounting. It is tested by
  \c test/cmocka_scheduler.c.
- \c tools/scheduler_http.c &ndash; the built-in HTTP(S) and websocket
  server on the scheduler port.
- \c tools/scheduler_lws.c &ndash; the same with libwebsockets on port + 100.
- \c tools/scheduler_proxy.c &ndash; commands from the browser and the
  forwarding of \c maos data to it.
- \c sys/scheduler_client.c &ndash; the client side used by \c maos,
  \c monitor, and \c drawdaemon.

All sockets are served by one thread with \c poll() in \c listen_port()
(\c sys/sock.c). A handler must therefore not block: replies are written with
a timeout and a client that does not read is disconnected.

\section sched-conn Connections

The scheduler listens on a unix socket \c $TEMP/scheduler and on TCP port
\c PORT (\c ~/.aos/port, or computed from the user name). \c maos on the same
host uses the unix socket. Connections from other hosts, \c monitor and
\c drawdaemon connecting by host name, and \c maos_client.py use TCP.

The unix socket is protected by file permissions and is not authenticated.
A TCP connection must first authenticate with \c CMD_AUTH; any other command
on an unauthenticated TCP connection closes it without a reply.

\c CMD_AUTH is a challenge-response: the scheduler sends a random nonce
string, the client replies with the string \c base64(sha1(nonce+password)),
and the scheduler replies with int 0 (authenticated for the lifetime of the
connection) or -1 (and closes). The password is the first line of
\c ~/.aos/htpasswd. It protects against replaying a recorded answer, but the
traffic is not encrypted and a man in the middle can take over an
authenticated connection. Use the TCP port only on a trusted network, or
tunnel it over SSH.

\section sched-format Message Format

All integers are native (little endian on supported platforms) 32-bit
\c int. A string is sent as an int length, including the terminating NUL,
followed by the bytes (\c stwritestr() and \c streadstr(); strings longer than
1 MiB are rejected). A file descriptor is passed with \c SCM_RIGHTS
(\c stwritefd()) and works over the unix socket only.

Every request starts with \c int[2] = {command, arg}. The meaning of \c arg
depends on the command; for commands sent by \c maos about its own job it is
the pid of \c maos and must be positive. The table lists the payload that
follows the header and the reply. Commands without a reply get none.

| Command | Value | Sent by | arg | Payload | Reply |
|---|---|---|---|---|---|
| CMD_START   | 1  | maos | pid | int nthread, int flag: bit 0 to wait for the scheduler, bits 1 and up the number of GPUs | int S_START when it may start, if waiting |
| CMD_FINISH  | 2  | maos | pid | | |
| CMD_STATUS  | 3  | maos | pid | status_t | |
| CMD_CRASH   | 4  | maos | pid | | |
| CMD_MONITOR | 5  | monitor | option: INT_MIN with bit 0 for load, bit 1 for plotting, bit 2 for text | | the job list, then updates (see below) |
| CMD_PATH    | 6  | maos | pid | string: command line with \c \\n between arguments | |
| CMD_KILL    | 7  | monitor | job pid | | |
| CMD_TRACE   | 8  | maos | pid | string: \c addr2line command | string |
| CMD_PROBE   | 9  | monitor | 0 | | int[3] {MON_VERSION, 0, 0} |
| CMD_DRAWCLI | 10 | draw() | >0: save the fd with this id; -10000..-1: get the fd saved with id -arg; otherwise request a new drawdaemon for job arg+10000 | fd if saving | int 0 and fd, or int -1 |
| CMD_REMOVE  | 11 | monitor | job pid, or -1 for all | | |
| CMD_DRAWSER | 12 | drawdaemon, monitor | >0: job pid; <=0: -id of a pending draw() request | | int 0 or -1 |
| CMD_MAOSCLI | 13 | maos_client.py | job pid | | int 0 or -1, then the connection is passed to maos |
| CMD_MAOSSER | 14 | maos | pid | | (keeps the connection to send MAOS_* commands) |
| CMD_RESTART | 15 | monitor | job pid | | |
| CMD_KILLED  | 16 | maos | pid | | |
| CMD_DUMMY   | 17 | (unused) | 0 | | |
| CMD_LAUNCH  | 18 | maos -d, skyc -d | number of strings: 2 for exe and command line, 1 for the command line only | strings | int 0 or -1 |
| CMD_AUTH    | 19 | TCP clients | 0 | see above | |

CMD_DRAWCLI passes descriptors and is only accepted over the unix socket.
CMD_DRAWSER and CMD_MAOSCLI hand the connection to another process and
remove it from the scheduler.

\section sched-monitor Messages to Monitors

The GTK monitor receives \c int[3] = {MON_*, pidnew, pid} messages, with
\c pidnew different from \c pid when a queued job is launched or restarted:

- \c MON_STATUS (3) followed by \c status_t.
- \c MON_PATH (6) followed by the job path as a string.
- \c MON_LOAD (14) with the third int holding the memory usage percentage in
  the lower 16 bits and the CPU usage percentage in the upper bits.
- \c MON_DRAWDAEMON (2) with the job pid asks the monitor to start a
  drawdaemon.
- \c MON_VERSION (13) answers CMD_PROBE.

The browser (websocket protocol \c maos-monitor-protocol) receives text
messages instead, each terminated by \c $ and with fields separated by \c &:
- \c pid&PATH&path$
- \c pid&STATUS&pidnew&info&start time&errhi&errlo&iseed&nseed&isim&nsim&rest&total&step time$
- \c pid&DRAW$ to open a plot.

Each command it sends is \c pid&NAME followed by \c ; (or a newline), where
NAME is \c REMOVE, \c KILL, \c RESTART, \c DRAW, \c DRAW_FIGFN, \c DRAW_PAUSE, \c DRAW_RESUME, or
\c DRAW_SINGLE (see \c ws_proxy_command()). Plot data from \c maos is
forwarded to the browser as binary messages through the
\c maos-drawdaemon-protocol connection.

\section sched-test Testing

\c make \c check-scheduler in the \c test folder of the build directory
starts a scheduler on a temporary home folder and port and runs the Python
scripts in \c test/scheduler, which check the commands above over both
sockets, authentication, invalid input, and the websocket interfaces. It is
not part of \c make \c check because it starts processes and uses the
network.
