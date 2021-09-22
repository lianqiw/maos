#!env python3
import aolib
import readbin
import socket
import sys
"""Get variables from running MAOS"""

sock=None #saves the socket
host=None
port=None
pid=0
def _pack_int(var):
    return var.to_bytes(4, byteorder='little')
def _pack_str(var):
    if len(var)>0:
        return _pack_int(len(var))+var.encode()
    else:
        return b''

def _recv_int(s):
    return int.from_bytes(s.recv(4), byteorder='little', signed=True)

def pause(pause):
    """Pause or unpause maos execution"""
    if not sock:
        print('Please run connect(host, port, [pid]) first')
        return
    msg=_pack_int(21)+_pack_int(pause)
    sock.send(msg)

def get_var(name):
    """Get value of an variable from maos with name"""
    if not sock:
        print('Please run connect(host, port, [pid]) first')
        return None
    msg=_pack_int(20)+_pack_int(1) 
    ans=sock.send(msg)
    msg=_pack_str(name)
    ans=sock.send(msg)
    ans=readbin.readbin(sock) #ans=aolib.readsock(sock.fileno())
    print(name, ans.shape, ans.dtype)
    return ans

def get_list():
    """Get the list of all variables"""
    if not sock:
        print('Please run connect(host, port, [pid]) first')
        return None
    msg=_pack_int(20)+_pack_int(1) 
    ans=sock.send(msg)
    msg=_pack_str("list")
    ans=sock.send(msg)
    (ans, header)=readbin.readbin(sock, 1) #ans=aolib.readsock(s.fileno())
    ans=header['COMMENT'].rstrip("\n\0").split("\n")[1:]
    return ans

def get_all():
    """Get value of all variables from maos"""
    vars=get_list()
    res={}
    for nf in vars:
        print(nf)
        res[nf]=get_var(nf)

    return res

def connect(host_=None, port_=None, pid_=None):
    """Connect to server that is running maos"""
    global sock, host, port, pid
    if host_ is not None:
        host=host_
    if port_ is not None:
        port=port_
    elif port is None:
        port=10000 #default port
    if pid_ is not None:
        pid=pid_

    if host is None:
        raise(Exception('Please call connect with host first'))
    if sock is not None:
        try:
            sock.close()
            sock=None
        except:
            pass

    try:
        sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
        sock.connect((host, port))
        msg=_pack_int(13)+_pack_int(pid)
        sock.send(msg)
        ans=_recv_int(sock)
        if ans==-1:
            sock.close()
            print('connect to maos failed')
            sock=None
    except Exception as error:
        print('connect/send failed:', error)
    
    
