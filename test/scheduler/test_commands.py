#!/usr/bin/env python3
"""Commands with invalid input must not crash or hang the scheduler. Usage: test_commands.py HOME"""
from common import *

# 1. CMD_DRAWCLI passes descriptors and is only valid over the unix socket.
s = tcp(); s.sendall(ints(CMD_DRAWCLI, 1) + b"x"); check(closed(s), "drawcli over tcp closes"); s.close()
s = unix(); s.sendall(ints(CMD_DRAWCLI, 1) + b"x"); check(closed(s), "drawcli without fd closes"); s.close()
s = unix(); s.sendall(ints(CMD_DRAWCLI, -1)); check(recv_int(s) == -1, "drawcli without saved drawdaemon answers -1"); s.close()
alive("drawcli")

# 2. CMD_TRACE runs the command with popen(), so only addr2line without shell characters is allowed.
flag = HOME + "/pwned"
for cmd in ["addr2line -f -e /bin/true; touch " + flag, "touch " + flag, "addr2line -f -e /bin/true 0x1$(touch " + flag + ")"]:
    s = unix(); s.sendall(ints(CMD_TRACE, os.getpid()) + sstr(cmd)); closed(s); s.close()
check(not os.path.exists(flag), "trace injection rejected")
s = unix(); s.sendall(ints(CMD_TRACE, os.getpid()) + sstr("addr2line -f -e /bin/true 0x1000 0x2000"))
n = recv_int(s); ans = recvn(s, n); s.close()
check(n > 0, "valid trace answered: %r" % ans)
s = unix(); s.sendall(ints(CMD_TRACE, os.getpid()) + ints(0)); check(closed(s), "trace with empty string"); s.close()
alive("trace")

# 3. CMD_LAUNCH edge cases
s = tcp(); s.sendall(ints(CMD_LAUNCH, 2) + sstr("/x/maos") + ints(0)); r = recv_int(s); s.close()
check(r == -1, "launch with empty command rejected")
s = tcp(); s.sendall(ints(CMD_LAUNCH, 1) + ints(0)); r = recv_int(s); s.close()
check(r == -1, "old launch with empty command rejected")
s = tcp(); s.sendall(ints(CMD_LAUNCH, 2) + ints(0x7fffffff)); r = recv_int(s)
check(r == -1 and closed(s), "launch with huge string length rejected"); s.close()
alive("launch")

# 4. commands from maos about its own job need a positive pid
for c in (CMD_START, CMD_PATH, CMD_STATUS, CMD_FINISH, CMD_MAOSSER):
    s = unix(); s.sendall(ints(c, -1, 1, 1)); check(closed(s), "cmd %d with pid -1 rejected" % c); s.close()
alive("negative pid")

# 5. CMD_PATH with a failing string read must not leave a dangling path.
mon = tcp(); mon.sendall(ints(CMD_MONITOR, -1))  # load reporting
time.sleep(0.2)
s = unix(); s.sendall(ints(CMD_PATH, os.getpid()) + sstr("first path")); time.sleep(0.2)
s.sendall(ints(CMD_PATH, os.getpid()) + ints(100) + b"short"); s.close(); time.sleep(0.3)
s = unix(); s.sendall(ints(CMD_PATH, os.getpid()) + sstr("second path")); time.sleep(0.3)
s.sendall(ints(CMD_FINISH, os.getpid())); time.sleep(0.3); s.close()
alive("path failure")
# a new monitor gets the job list, including the finished job
m2 = tcp(); m2.sendall(ints(CMD_MONITOR, 53)); m2.settimeout(2)
seen = {}
try:
    while True:
        h = struct.unpack("<3i", recvn(m2, 12))
        if h[0] == MON_PATH:
            seen.setdefault(h[2], []).append(recv_str(m2).rstrip(b"\0").decode())
        elif h[0] == MON_STATUS:
            recvn(m2, STATUS_SIZE)
except Exception:
    pass
m2.close(); mon.close()
check(seen.get(os.getpid()) == ["second path"], "job path after failed update: %r" % seen.get(os.getpid()))
s = tcp(); s.sendall(ints(CMD_REMOVE, os.getpid())); s.close()
alive("monitor")

# 6. A job command socket is closed and its descriptor reused by another
# client. maos_command() must not write to the new connection.
s = unix(); s.sendall(ints(CMD_START, os.getpid(), 1, 0)); time.sleep(0.2)
c = unix(); c.sendall(ints(CMD_MAOSSER, os.getpid())); time.sleep(0.2); c.close(); time.sleep(0.3)
m = tcp(); probe(m)  # likely reuses the descriptor of c
cli = tcp(); cli.sendall(ints(CMD_MAOSCLI, os.getpid())); r = cli.recv(4); cli.close()
check(r == ints(-1), "maos_command after sock_cmd closed fails: %r" % r)
m.settimeout(1)
try:
    extra = m.recv(100)
except socket.timeout:
    extra = b""
check(extra == b"", "no stray command on reused socket: %r" % extra)
m.close()
s.sendall(ints(CMD_FINISH, os.getpid())); s.close(); time.sleep(0.2)
s = tcp(); s.sendall(ints(CMD_REMOVE, os.getpid())); s.close()
alive("stale sock")

# 7. invalid command and garbage
s = tcp(); s.sendall(ints(1000, 0)); check(closed(s), "invalid command closed"); s.close()
s = tcp(); s.sendall(b"\xff" * 64); check(closed(s), "garbage closed"); s.close()
s = tcp_raw(); s.sendall(b"GE"); s.close()
alive("garbage")
finish()
