#!/usr/bin/env python3
"""Get variables from a running MAOS through the scheduler.

Usage:
    import maos_client
    maos_client.connect('localhost', pid_=12345) #port defaults to the scheduler port
    names=maos_client.get_list()
    gradcl=maos_client.get_var('simu.gradcl')
"""
import os
import socket
import readbin

sock=None #saves the socket
fp=None   #buffered reader of sock
host=None
port=None
pid=0

CMD_MAOSCLI=13 #scheduler command to link a client to maos (sys/scheduler_client.h)
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
    data=b''
    while len(data)<4:
        chunk=s.recv(4-len(data))
        if not chunk:
            raise ConnectionError('connection closed')
        data+=chunk
    return int.from_bytes(data, byteorder='little', signed=True)

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
    if not sock:
        print('Please run connect(host, port, pid) first')
        return False
    return True

def _read_reply():
    '''Read a variable sent by maos (writesock) and the end of data marker that follows it.'''
    (out, err)=readbin.readbin_auto(fp, False)
    magic=readbin.readbin_magic(fp)
    if magic!=readbin.M_EOD:
        print('maos_client: expected end of data, got', hex(magic))
    if err<0:
        print('maos_client: read failed', err)
    return out

def pause(pause):
    """Pause (1) or resume (0) maos execution"""
    if not _check():
        return
    sock.sendall(_pack_int(MAOS_PAUSE)+_pack_int(pause))

def get_var(name):
    """Get value of an variable from maos with name, such as simu.gradcl"""
    if not _check():
        return None
    sock.sendall(_pack_int(MAOS_VAR)+_pack_int(1)+_pack_str(name))
    return _read_reply()

def get_list():
    """Get the list of all variables"""
    if not _check():
        return None
    sock.sendall(_pack_int(MAOS_VAR)+_pack_int(1)+_pack_str("list"))
    ans=_read_reply()
    header=getattr(ans, 'header', '')
    if isinstance(header, str):
        header=header.split('\n')
    #The first line is "Available variables are:"
    return [str(line).strip() for line in header[1:] if str(line).strip()]

def get_all():
    """Get value of all variables from maos"""
    vars=get_list()
    res={}
    if vars is None:
        return res
    for nf in vars:
        res[nf]=get_var(nf)
    return res

def connect(host_=None, port_=None, pid_=None):
    """Connect to the scheduler on host that runs maos with pid. The port defaults to the scheduler port."""
    global sock, fp, host, port, pid
    if host_ is not None:
        host=host_
    if port_ is not None:
        port=port_
    elif port is None:
        port=default_port()
    if pid_ is not None:
        pid=pid_

    if host is None:
        raise(Exception('Please call connect with host first'))
    if sock is not None:
        try:
            fp.close()
            sock.close()
        except OSError:
            pass
        sock=None
        fp=None

    try:
        sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((host, port))
        sock.sendall(_pack_int(CMD_MAOSCLI)+_pack_int(pid))
        ans=_recv_int(sock)
        if ans==-1:
            sock.close()
            print('connect to maos failed')
            sock=None
        else:
            fp=sock.makefile('rb')
    except Exception as error:
        print('connect/send failed:', error)
        sock=None
        fp=None
