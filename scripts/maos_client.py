#!env python3
import aolib
import readbin
#get variables from running MAOS
import socket
s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)    
s.connect(('maxwell', 10000))

def send_int(var):
    return var.to_bytes(4, byteorder='little')
def send_str(var):
    if len(var)>0:
        return send_int(len(var))+var.encode()
    else:
        return b''



def recv_int():
    return int.from_bytes(s.recv(4), byteorder='little')



def maos_pause(pause):
    msg=send_int(21)+send_int(pause)
    s.send(msg)

def maos_get_var(name):
    msg=send_int(20)+send_int(1) 
    ans=s.send(msg)
    msg=send_str(name);
    ans=s.send(msg)
    ans=readbin.readbin(s) #ans=aolib.readsock(s.fileno())
    print(name, ans.shape, ans.dtype)
    return ans

def maos_get_list():
    msg=send_int(20)+send_int(1) 
    ans=s.send(msg)
    msg=send_str("list");
    ans=s.send(msg)
    (ans, header)=readbin.readbin(s, 1) #ans=aolib.readsock(s.fileno())
    ans=header.rstrip("\n\0").split("\n")[1:]
    return ans


pid=0
msg=send_int(13)+send_int(pid)
s.send(msg)
ans=recv_int()

vars=maos_get_list()
res={}
for nf in vars :
    res[nf]=maos_get_var(nf)

s.close()

