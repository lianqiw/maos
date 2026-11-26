#!/usr/bin/env python3
"""Get variables from a running MAOS through the scheduler.

The client connects to the scheduler over TCP, authenticates with the password
in ~/.aos/htpasswd (which must match the one on the scheduler host), and asks
the scheduler (CMD_MAOSCLI) to hand the connection to the maos process with the
given pid. maos then serves requests on that connection in maos_var()
(maos/maos.c):

- MAOS_VAR, 1, name: maos replies with the variable in .bin format (writesock)
  followed by an M_EOD marker. An unknown name (or "list") returns an empty
  cell whose header lists the available names (e.g. "simu.gradcl",
  "powfs[0].saloc").
- MAOS_PAUSE, flag: pause (1) or resume (0) the simulation. No reply.

Only one connection is kept, in module level variables. Requires readbin.py
(same directory) and numpy/scipy.

Usage:
    import maos_client
    maos_client.connect('localhost', pid_=12345) #port defaults to the scheduler port
    names=maos_client.get_list()
    gradcl=maos_client.get_var('simu.gradcl')
    maos_client.pause(1)  #pause; pause(0) resumes
    maos_client.close()
"""
import os
import socket
import hashlib
import base64
import readbin

sock=None #saves the socket
fp=None   #buffered reader of sock
host=None
port=None
pid=0
timeout=None #socket timeout in seconds used by the current connection

CMD_MAOSCLI=13 #scheduler command to link a client to maos (sys/scheduler_client.h)
CMD_AUTH=19    #scheduler command to authenticate a TCP connection
MAOS_VAR=20    #maos command to get or put a variable
MAOS_PAUSE=21  #maos command to pause or resume

def _pack_int(var):
    return int(var).to_bytes(4, byteorder='little', signed=True)
def _pack_str(var):
    '''Pack a string as stwritestr() does: length including the terminating null, then the string.'''
    if var:
        data=var.encode()+b'\0'
        return _pack_int(len(data))+data
    else:
        return _pack_int(0)

def _recv_int(s):
    '''Receive a 4 byte little endian signed integer. Raises ConnectionError if the peer closes.'''
    data=b''
    while len(data)<4:
        chunk=s.recv(4-len(data))
        if not chunk:
            raise ConnectionError('connection closed')
        data+=chunk
    return int.from_bytes(data, byteorder='little', signed=True)

def _recv_str(s):
    '''Receive a string written by stwritestr().'''
    n=_recv_int(s)
    data=b''
    while len(data)<n:
        chunk=s.recv(n-len(data))
        if not chunk:
            raise ConnectionError('connection closed')
        data+=chunk
    return data.rstrip(b'\0')

def _auth(s):
    '''Authenticate with the scheduler as scheduler_auth() does, using the password in ~/.aos/htpasswd.'''
    fn=os.path.expanduser('~/.aos/htpasswd')
    try:
        with open(fn, 'rb') as f:
            passwd=f.readline().strip()
    except OSError as error:
        raise ConnectionError(f'{fn} is not readable ({error}). Copy it from the host that runs the scheduler.') from error
    s.sendall(_pack_int(CMD_AUTH)+_pack_int(0))
    nonce=_recv_str(s)
    key=base64.b64encode(hashlib.sha1(nonce+passwd).digest()).decode()
    s.sendall(_pack_str(key))
    if _recv_int(s)!=0:
        raise ConnectionError('authentication failed: ~/.aos/htpasswd must be the same as on the scheduler host')

def _hashlittle(key, initval=0):
    '''Port of hashlittle() in sys/hashlittle.c for byte strings.'''
    M=0xFFFFFFFF
    def rot(x, k):
        return ((x<<k)|(x>>(32-k)))&M
    length=len(key)
    dstep=max(1, length//10000)
    a=b=c=(0xdeadbeef+length+initval)&M
    data=key+b'\0'*12 #masked reads beyond the end see zeros
    def word(i):
        return int.from_bytes(data[4*i:4*i+4], 'little')
    k=0
    while length>dstep*12:
        a=(a+word(k))&M; b=(b+word(k+1))&M; c=(c+word(k+2))&M
        a=(a-c)&M; a^=rot(c, 4);  c=(c+b)&M
        b=(b-a)&M; b^=rot(a, 6);  a=(a+c)&M
        c=(c-b)&M; c^=rot(b, 8);  b=(b+a)&M
        a=(a-c)&M; a^=rot(c, 16); c=(c+b)&M
        b=(b-a)&M; b^=rot(a, 19); a=(a+c)&M
        c=(c-b)&M; c^=rot(b, 4);  b=(b+a)&M
        length-=dstep*12
        k+=dstep*3
    if length==0:
        return c
    def part(i, n):#first n bytes of word i
        return int.from_bytes(data[4*i:4*i+n], 'little')
    if length<=12: #like the C switch, a longer remainder from sampling is skipped
        a=(a+part(k, min(length, 4)))&M
        if length>4:
            b=(b+part(k+1, min(length-4, 4)))&M
        if length>8:
            c=(c+part(k+2, length-8))&M
    c^=b; c=(c-rot(b, 14))&M
    a^=c; a=(a-rot(c, 11))&M
    b^=a; b=(b-rot(a, 25))&M
    c^=b; c=(c-rot(b, 16))&M
    a^=c; a=(a-rot(c, 4))&M
    b^=a; b=(b-rot(a, 14))&M
    c^=b; c=(c-rot(b, 24))&M
    return c

def default_port():
    '''Scheduler port as chosen by init_hosts() in sys/scheduler_client.c: ~/.aos/port, or a hash of the user name.'''
    try:
        with open(os.path.expanduser('~/.aos/port')) as f:
            p=int(f.read().split()[0])
            if p>=1000:
                return p
    except (OSError, ValueError, IndexError):
        pass
    user=os.environ.get('USER') or os.environ.get('USERNAME') or ''
    return (_hashlittle(user.encode())&0x2FFF)|10000

def _check():
    '''Return True if connected, otherwise print a hint and return False.'''
    if not sock:
        print('Please run connect(host, port, pid) first')
        return False
    return True

def _read_reply():
    '''Read a variable sent by maos (writesock) and the end of data marker that follows it.'''
    (out, err)=readbin.readbin_auto(fp, False)
    if err==1: #M_EOD read in place of the data: maos sent nothing
        return None
    magic=readbin.readbin_magic(fp)
    if magic!=readbin.M_EOD:
        print('maos_client: expected end of data, got', hex(magic), '; the connection may be out of sync, consider reconnecting.')
    if err<0:
        print('maos_client: read failed', err)
    return out

def _close():
    '''Close the connection, ignoring errors. The host, port, and pid are kept for connect().'''
    global sock, fp
    for obj in (fp, sock):
        if obj is not None:
            try:
                obj.close()
            except OSError:
                pass
    sock=None
    fp=None

def close():
    """Close the connection to maos. connect() with no arguments reconnects."""
    _close()

def _request(data, reply=True):
    '''Send a request and optionally read the reply. On socket errors (including
    timeout) the connection is closed, since the stream is no longer in sync,
    and the error is raised.'''
    try:
        sock.sendall(data)
        return _read_reply() if reply else None
    except OSError:
        _close()
        raise

def pause(pause):
    """Pause (1) or resume (0) maos execution. maos does not reply.

    Raises OSError (and closes the connection) if sending fails."""
    if not _check():
        return
    _request(_pack_int(MAOS_PAUSE)+_pack_int(pause), reply=False)

def _header_lines(ans):
    '''Non-empty lines of the header of ans. readbin gives the header as a str, or an array of lines if multi-line.'''
    header=getattr(ans, 'header', '')
    if isinstance(header, str):
        header=header.split('\n')
    return [str(line).strip() for line in header if str(line).strip()]

def get_var(name):
    """Get the value of a variable from maos by name, such as 'simu.gradcl' or 'powfs[0].saloc'.

    Returns a numpy array (or object array for cells, or scipy sparse matrix)
    as returned by readbin, with the .bin header in its .header attribute, or
    None if not connected or maos sent nothing. For an unknown name maos
    returns an empty array whose header lists the available variables; a
    warning is printed in that case. Raises OSError (and closes the
    connection) on socket errors or timeout. Use get_list() for valid names."""
    if not _check():
        return None
    if not name:
        raise ValueError('get_var: name must be a non-empty string')
    out=_request(_pack_int(MAOS_VAR)+_pack_int(1)+_pack_str(name))
    if name!='list' and _header_lines(out)[:1]==['Available variables are:']:
        print(f'maos_client: {name} is not found in maos. Use get_list() for valid names.')
    return out

def get_list():
    """Get the list of all variable names that get_var() accepts, or None if not connected."""
    if not _check():
        return None
    ans=_request(_pack_int(MAOS_VAR)+_pack_int(1)+_pack_str("list"))
    #The first line is "Available variables are:"
    return _header_lines(ans)[1:]

def get_all():
    """Get the value of all variables from maos as a dict keyed by name. Empty if not connected."""
    vars=get_list()
    res={}
    if vars is None:
        return res
    for nf in vars:
        res[nf]=get_var(nf)
    return res

def connect(host_=None, port_=None, pid_=None, timeout_=30):
    """Connect to the scheduler on host that runs maos with pid and link to that maos.

    host_:    scheduler host name or address (TCP only; required on first call).
    port_:    scheduler port. Defaults to default_port().
    pid_:     pid of the running maos.
    timeout_: socket timeout in seconds for connecting and for each subsequent
              send/receive (not the whole transfer). None blocks forever.

    Arguments that are None reuse the values of the previous call, so connect()
    reconnects. Any existing connection is closed first. On failure, a message
    is printed and the module stays disconnected (no exception is raised,
    except when no host is known)."""
    global sock, fp, host, port, pid, timeout
    if host_ is not None:
        host=host_
    if port_ is not None:
        port=port_
    elif port is None:
        port=default_port()
    if pid_ is not None:
        pid=pid_

    timeout=timeout_

    if host is None:
        raise(Exception('Please call connect with host first'))
    _close()

    try:
        #create_connection resolves the host and applies the timeout to connect and later operations.
        sock=socket.create_connection((host, port), timeout=timeout)
        _auth(sock)
        sock.sendall(_pack_int(CMD_MAOSCLI)+_pack_int(pid))
        ans=_recv_int(sock) #scheduler replies 0 if the connection is passed to maos, -1 otherwise.
        if ans!=0:
            _close()
            print(f'connect to maos failed: scheduler at {host}:{port} cannot reach maos with pid {pid}')
        else:
            fp=sock.makefile('rb')
    except (OSError, ValueError) as error: #ConnectionError and socket.timeout are OSError
        print(f'connect/send to {host}:{port} failed:', error)
        _close()
