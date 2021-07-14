#!env python3
import aolib
import readbin
import socket
import sys
"""Get variables from running MAOS"""

s=None #saves the socket
host=None
port=None
def _pack_int(var):
    return var.to_bytes(4, byteorder='little')
def _pack_str(var):
    if len(var)>0:
        return _pack_int(len(var))+var.encode()
    else:
        return b''

def _recv_int(s):
    return int.from_bytes(s.recv(4), byteorder='little')

def pause(pause):
    """Pause or unpause maos execution"""
    if not s:
        print('Please run connect(host, port, [pid]) first')
        return None
    msg=_pack_int(21)+_pack_int(pause)
    s.send(msg)

def get_var(name):
    """Get value of an variable from maos with name"""
    if not s:
        print('Please run connect(host, port, [pid]) first')
        return None
    msg=_pack_int(20)+_pack_int(1) 
    ans=s.send(msg)
    msg=_pack_str(name)
    ans=s.send(msg)
    ans=readbin.readbin(s) #ans=aolib.readsock(s.fileno())
    print(name, ans.shape, ans.dtype)
    return ans

def get_list():
    """Get the list of all variables"""
    if not s:
        print('Please run connect(host, port, [pid]) first')
        return None
    msg=_pack_int(20)+_pack_int(1) 
    ans=s.send(msg)
    msg=_pack_str("list")
    ans=s.send(msg)
    (ans, header)=readbin.readbin(s, 1) #ans=aolib.readsock(s.fileno())
    ans=header['COMMENT'].rstrip("\n\0").split("\n")[1:]
    return ans

def get_all():
    """Get value of all variables from maos"""
    if not s:
        print('Please run connect(host, port, [pid]) first')
        return None
    vars=get_list()
    res={}
    for nf in vars:
        res[nf]=get_var(nf)

    return res

def connect(host_, port_, pid=0):
    """Connect to server that is running maos"""
    global s, host, port
    host=host_
    port=port_
    
    s=socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.connect((host, port))
    
    msg=_pack_int(13)+_pack_int(pid)
    s.send(msg)
    ans=_recv_int(s)
    #pause(1)
    
