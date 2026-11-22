#!/usr/bin/env python3
"""TCP clients must authenticate with CMD_AUTH. Usage: test_auth.py HOME BINDIR PROBE SRCDIR
PROBE is the scheduler_probe program, which uses scheduler_connect()."""
import shutil, subprocess
from common import *
BIN, PROBE, SRCDIR = sys.argv[2], sys.argv[3], sys.argv[4]

# 1. commands over TCP without authentication are rejected without any action or reply
for c in range(1, 19):
    s = tcp_raw(); s.sendall(ints(c, os.getpid()) + b"\0" * 16)
    r = reply(s); s.close()
    check(r == b"", "unauthenticated cmd %d is rejected (%r)" % (c, r))
s = tcp_raw(); s.sendall(ints(CMD_LAUNCH, 2) + sstr("/bin/sh") + sstr("/tmp/x/maos -c true x"))
check(reply(s) == b"", "unauthenticated LAUNCH gets no answer"); s.close()
# 2. authentication
s = tcp_raw(); check(auth(s) == 0 and probe(s), "authenticated PROBE works")
check(probe(s), "authentication lasts for the connection"); s.close()
s = tcp_raw(); r = auth(s, passwd=b"wrong"); check(r == -1 and reply(s) == b"", "wrong password is rejected and closed"); s.close()
s = tcp_raw(); check(auth(s, passwd=b"") == -1, "empty password is rejected"); s.close()
s = tcp_raw(); s.sendall(ints(CMD_AUTH, 0)); recv_str(s)
key = b"AAAAAAAAAAAAAAAAAAAAAAAAAAA=\0"; s.sendall(ints(len(key)) + key)
check(recv_int(s) == -1, "wrong answer is rejected"); s.close()
s = tcp_raw(); s.sendall(ints(CMD_AUTH, 0)); s.close()  # client disconnects in the middle
s = tcp_raw(); s.sendall(ints(CMD_AUTH, 0)); recv_str(s)
s.sendall(ints(1 << 30)); check(recv_int(s) == -1, "huge answer is rejected"); s.close()
# 3. the unix socket does not need authentication
s = unix(); check(probe(s), "unix socket PROBE without authentication"); s.close()
# 4. the C client: scheduler_connect() authenticates
env = dict(os.environ, HOME=HOME)
r = subprocess.run([PROBE, "localhost"], env=env, capture_output=True, text=True, timeout=30)
check("probe ok 13" in r.stdout, "C client authenticates: %s" % r.stdout.strip())
other = HOME + "/other"; shutil.rmtree(other, ignore_errors=True); os.makedirs(other + "/.aos")
shutil.copy(HOME + "/.aos/port", other + "/.aos/port")
env_other = dict(os.environ, HOME=other)
r = subprocess.run([PROBE, "localhost"], env=env_other, capture_output=True, text=True, timeout=30)
check("connect failed" in r.stdout and "htpasswd is not found" in r.stderr + r.stdout, "C client without htpasswd fails with a message")
open(other + "/.aos/htpasswd", "w").write("wrongpassword\n")
r = subprocess.run([PROBE, "localhost"], env=env_other, capture_output=True, text=True, timeout=30)
check("connect failed" in r.stdout and "Authentication failed" in r.stderr + r.stdout, "C client with a wrong password fails with a message")
# 5. maos_client.py authenticates. The maos pid does not exist, so the scheduler answers -1 after authentication.
code = ("import sys, types; sys.modules['readbin'] = types.ModuleType('readbin'); sys.path.insert(0, %r); "
        "import maos_client; maos_client.connect('localhost', %d, 1)" % (SRCDIR + "/scripts", PORT))
r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=30)
check("connect to maos failed" in r.stdout and "authentication" not in r.stdout, "maos_client.py authenticates: %r" % r.stdout.strip())
r = subprocess.run([sys.executable, "-c", code], env=env_other, capture_output=True, text=True, timeout=30)
check("authentication failed" in r.stdout, "maos_client.py with a wrong password: %r" % r.stdout.strip())
# 6. a detached maos is queued through the scheduler over TCP (scheduler_launch_exe)
if os.path.exists(BIN + "/maos"):
    out = HOME + "/detached"
    r = subprocess.run([BIN + "/maos", "-d", "-o", out, "aper.d=5", "sim.end=10", "-g-1", "-n", "1"], env=env, cwd=HOME,
                       capture_output=True, text=True, timeout=60)
    time.sleep(1)
    m = tcp(); m.sendall(ints(CMD_MONITOR, 53)); m.settimeout(2)
    paths = []
    try:
        while True:
            h = struct.unpack("<3i", recvn(m, 12))
            if h[0] == MON_PATH:
                paths.append(recv_str(m).decode(errors="replace"))
            elif h[0] == MON_STATUS:
                recvn(m, STATUS_SIZE)
    except Exception:
        pass
    m.close()
    check(any("detached" in p for p in paths), "detached maos is queued through the scheduler (rc=%d)" % r.returncode)
s = tcp_raw(); s.sendall(b"\xff" * 64); check(reply(s) == b"", "garbage closed"); s.close()
s = unix(); check(probe(s), "alive"); s.close()
finish()
