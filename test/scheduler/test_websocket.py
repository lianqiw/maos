#!/usr/bin/env python3
"""Test the browser interface over websocket. Usage: test_websocket.py HOME [lws]
Without lws, the built-in server on PORT is tested. With lws, libwebsockets on PORT+100."""
import subprocess
from common import *
LWS = len(sys.argv) > 2 and sys.argv[2] == "lws"
WSPORT = PORT + 100 if LWS else PORT
closed_ws = False

def ws_read(s):
    h = recvn(s, 2); n = h[1] & 0x7f
    if n == 126: n = struct.unpack(">H", recvn(s, 2))[0]
    elif n == 127: n = struct.unpack(">Q", recvn(s, 8))[0]
    return h[0], recvn(s, n)
def ws_msgs(s, t):
    """Read complete messages (reassembling fragments) for up to t seconds."""
    global closed_ws
    out = []; cur = b""; end = time.time() + t
    s.settimeout(0.3)
    while time.time() < end:
        try:
            h = recvn(s, 2)
        except socket.timeout:
            continue
        except (EOFError, ConnectionResetError):
            closed_ws = True
            break
        try:
            n = h[1] & 0x7f
            if n == 126: n = struct.unpack(">H", recvn(s, 2))[0]
            elif n == 127: n = struct.unpack(">Q", recvn(s, 8))[0]
            s.settimeout(5); p = recvn(s, n); s.settimeout(0.3)
        except (EOFError, ConnectionResetError):
            closed_ws = True
            break
        cur += p
        if h[0] & 0x80:
            out.append((h[0] & 0x0f, cur)); cur = b""
    return out
def statuses(msgs, pid):
    """isim and info of the status messages of job pid."""
    isims = []; infos = []
    for op, p in msgs:
        for e in p.decode(errors="replace").split("$"):
            f = e.split("&")
            if len(f) > 10 and f[0] == str(pid) and f[1] == "STATUS":
                infos.append(int(f[3])); isims.append(int(f[9]))
    return isims, infos
def remove(pid):
    s = tcp(); s.sendall(ints(CMD_REMOVE, pid)); s.close()

if not LWS:
    # 0. authentication and frame parsing of the built-in server
    if not os.path.exists(HOME + "/.aos/htpasswd"):
        s, r = ws_connect(authorize=False); s.close()  # creates the password file
    s, r = ws_connect(authorize=False)
    check(r.startswith(b"HTTP/1.1 401"), "websocket without authorization gets 401"); s.close()
    s, r = ws_connect()
    check(r.startswith(b"HTTP/1.1 101"), "websocket upgrade")
    # a ping and a fragmented message in one write are both handled
    s.sendall(frame(b"hi", opcode=9) + frame(b"0&RES", fin=False) + frame(b"TART;", opcode=0) + frame(b"", opcode=9))
    pongs = []
    s.settimeout(2)
    try:
        while len(pongs) < 2:
            h, p = ws_read(s)
            if h == 0x8A: pongs.append(p)
    except Exception as e:
        print("     ws read:", e)
    check(pongs[:1] == [b"hi"] and len(pongs) == 2, "two pongs for pings around a fragmented message")
    s.sendall(frame(b"", length=1 << 40)[:14])
    check(closed(s), "oversized frame closes"); s.close()
    s, r = ws_connect(); s.sendall(frame(b"abc", mask=False)); check(closed(s), "unmasked frame closes"); s.close()
    s, r = ws_connect(); s.sendall(frame(b"abc", opcode=0)); check(closed(s), "unexpected continuation closes"); s.close()
    s = tcp_raw(); s.sendall(b"GET / HTTP/1.1\r\nUpgrade: websocket\r\n\r\n"); s.recv(200); s.close()
    alive("websocket frames")

# 1. the monitor receives the status of a job in order
job = subprocess.Popen(["sleep", "120"]); pid = job.pid
w, r = ws_connect(port=WSPORT)
check(r.startswith(b"HTTP/1.1 101"), "upgrade %r" % r[:20])
ws_msgs(w, 1)  # initial messages
m = unix(); m.sendall(ints(CMD_START, pid, 1, 0)); time.sleep(0.2)
m.sendall(b"".join(status(pid, i) for i in range(300)))
m.sendall(ints(CMD_FINISH, pid))
msgs = ws_msgs(w, 3)
isims, infos = statuses(msgs, pid)
check(set(op for op, p in msgs) == {1}, "all monitor messages are text frames")
check(isims == sorted(isims) and len(isims) > 0, "status in order (%d received)" % len(isims))
check(infos and infos[-1] == 11, "last status is finished: %r" % infos[-3:])
m.close(); job.kill(); remove(pid)
alive("ordered status")

# 2. a browser that is slow to read. Either all messages arrive in order, or
# the scheduler closes it and a new connection receives the full status.
job2 = subprocess.Popen(["sleep", "120"]); pid2 = job2.pid
w2, r = ws_connect(port=WSPORT)
w2.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 2048)
ws_msgs(w2, 1)
m = unix(); m.sendall(ints(CMD_START, pid2, 1, 0)); time.sleep(0.2)
m.sendall(b"".join(status(pid2, i) for i in range(20000)))
m.sendall(ints(CMD_FINISH, pid2)); time.sleep(3)
closed_ws = False
msgs = []
for k in range(60):  # until the finished status arrives or the connection is closed
    msgs += ws_msgs(w2, 1)
    if closed_ws or any(b"&STATUS&%d&11&" % pid2 in p for op, p in msgs[-5:]): break
if closed_ws:
    print("     slow client was closed after %d messages; reconnecting" % len(msgs))
    w2.close(); w2, r = ws_connect(port=WSPORT); msgs = ws_msgs(w2, 2)
isims, infos = statuses(msgs, pid2)
check(isims == sorted(isims) and len(isims) > 0, "slow reader: status in order (%d received)" % len(isims))
check(infos and infos[-1] == 11, "slow reader: last status is finished: %r" % infos[-3:])
m.close(); job2.kill(); w2.close(); remove(pid2)
alive("slow reader")

# 3. a ping is answered and a command is handled
w.sendall(frame(b"hi", 9))
check(any(op == 10 and p == b"hi" for op, p in ws_msgs(w, 1)), "pong")
w.sendall(frame(b"0&RESTART;")); time.sleep(0.2)
alive("command")
w.close(); time.sleep(0.2)
alive("monitor close")

# 4. drawing request for a job that does not exist
d, r = ws_connect("maos-drawdaemon-protocol", port=WSPORT)
d.sendall(frame(b"%d&DRAW;" % pid)); time.sleep(0.5); d.close()
alive("draw for missing job")

# 5. a browser that stops reading while many statuses are sent
job = subprocess.Popen(["sleep", "120"]); pid = job.pid
w, r = ws_connect(port=WSPORT)
w.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 4096)
m = unix(); m.sendall(ints(CMD_START, pid, 1, 0))
t0 = time.time()
for k in range(30):
    m.sendall(b"".join(status(pid, i) for i in range(k * 1000, k * 1000 + 1000)))
alive("stalled client (%.1fs)" % (time.time() - t0))
m.sendall(ints(CMD_FINISH, pid)); m.close(); job.kill()
w.close(); time.sleep(0.5); remove(pid)
alive("stalled client closed")
finish()
