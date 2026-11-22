"""Helpers for the scheduler tests. See run.sh and docs/94_scheduler.md.

The scheduler under test runs with HOME set to a private folder that holds
.aos/port and .aos/htpasswd. Usage of each test: test_*.py HOME [args]
"""
import base64, glob, hashlib, os, socket, struct, sys, time

HOME = sys.argv[1]
PORT = int(open(HOME + "/.aos/port").read())
UNIX = glob.glob(HOME + "/.aos/tmp/*")[0] + "/scheduler"
(CMD_START, CMD_FINISH, CMD_STATUS, CMD_CRASH, CMD_MONITOR, CMD_PATH, CMD_KILL, CMD_TRACE, CMD_PROBE, CMD_DRAWCLI,
 CMD_REMOVE, CMD_DRAWSER, CMD_MAOSCLI, CMD_MAOSSER, CMD_RESTART, CMD_KILLED, CMD_DUMMY, CMD_LAUNCH, CMD_AUTH) = range(1, 20)
MON_STATUS, MON_PATH, MON_VERSION, MON_LOAD = 3, 6, 13, 14
STATUS_SIZE = 144  # sizeof(status_t) on 64 bit platforms
ok = True

def ints(*v): return struct.pack("<%di" % len(v), *v)
def sstr(s):
    b = s.encode() + b"\0"
    return ints(len(b)) + b
def recvn(s, n):
    b = b""
    while len(b) < n:
        c = s.recv(n - len(b))
        if not c: raise EOFError("closed after %d of %d" % (len(b), n))
        b += c
    return b
def recv_int(s): return struct.unpack("<i", recvn(s, 4))[0]
def recv_str(s): return recvn(s, recv_int(s))
def password(): return open(HOME + "/.aos/htpasswd", "rb").readline().strip()
def auth(s, passwd=None):
    """Authenticate as scheduler_auth() in sys/scheduler_client.c. Return the scheduler's answer (0 on success)."""
    s.sendall(ints(CMD_AUTH, 0))
    nonce = recv_str(s).rstrip(b"\0")
    if passwd is None:  # read after the nonce: the scheduler creates the file if missing
        passwd = password()
    key = base64.b64encode(hashlib.sha1(nonce + passwd).digest()) + b"\0"
    s.sendall(ints(len(key)) + key)
    return recv_int(s)
def tcp_raw(port=PORT): return socket.create_connection(("127.0.0.1", port), timeout=5)
def tcp():
    s = tcp_raw()
    assert auth(s) == 0
    return s
def unix():
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM); s.settimeout(5); s.connect(UNIX); return s
def check(cond, tag):
    global ok
    print(("PASS" if cond else "FAIL"), tag); ok &= bool(cond)
def closed(s, t=3):
    """Read until EOF and return True, or False if the connection stays open. Data before the close is ignored."""
    s.settimeout(t)
    try:
        while s.recv(4096):
            pass
        return True
    except ConnectionResetError:
        return True
    except socket.timeout:
        return False
def reply(s):
    """Data received before the connection is closed, or None if it stays open."""
    s.settimeout(10); b = b""
    try:
        while True:
            c = s.recv(4096)
            if not c: return b
            b += c
    except ConnectionResetError:
        return b
    except socket.timeout:
        return None
def probe(s):
    s.sendall(ints(CMD_PROBE, 0))
    return struct.unpack("<3i", recvn(s, 12))[0] == MON_VERSION
def alive(tag):
    try:
        s = tcp(); good = probe(s); s.close()
    except Exception as e:
        good = False
    check(good, "alive after " + tag)
def status(pid, isim, info=1):
    """A CMD_STATUS message."""
    st = bytearray(STATUS_SIZE)
    struct.pack_into("<iii", st, 72, 0, 1, 1)  # iseed nseed nthread
    struct.pack_into("<iii", st, 84, isim, 0, 100000)  # isim simstart simend
    struct.pack_into("<i", st, 120, info)
    return ints(CMD_STATUS, pid) + bytes(st)
def ws_connect(proto="maos-monitor-protocol", port=PORT, authorize=True):
    """Open a websocket. Return the socket and the http response header."""
    s = tcp_raw(port)
    key = base64.b64encode(os.urandom(16)).decode()
    hdr = ("GET / HTTP/1.1\r\nHost: x\r\nUpgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Key: %s\r\n"
           "Sec-WebSocket-Version: 13\r\nSec-WebSocket-Protocol: %s\r\n" % (key, proto))
    if authorize and os.path.exists(HOME + "/.aos/htpasswd"):
        user = os.environ.get("USER", "nobody").encode()
        hdr += "Authorization: Basic %s\r\n" % base64.b64encode(user + b":" + password()).decode()
    s.sendall((hdr + "\r\n").encode())
    resp = b""
    while b"\r\n\r\n" not in resp:
        c = s.recv(1)
        if not c: break
        resp += c
    return s, resp
def frame(payload, opcode=1, fin=True, mask=True, length=None):
    """A client websocket frame."""
    b = bytes([(0x80 if fin else 0) | opcode])
    n = len(payload) if length is None else length
    m = 0x80 if mask else 0
    if n < 126: b += bytes([m | n])
    elif n < 65536: b += bytes([m | 126]) + struct.pack(">H", n)
    else: b += bytes([m | 127]) + struct.pack(">Q", n)
    if mask:
        k = os.urandom(4); b += k
        payload = bytes(c ^ k[i % 4] for i, c in enumerate(payload))
    return b + payload
def finish():
    print("ALL PASS" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)
