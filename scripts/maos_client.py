#!env python3
import aolib
import readbin
#get variables from running MAOS
import socket
import sys
res={}

def send_int(var):
    return var.to_bytes(4, byteorder='little')
def send_str(var):
    if len(var)>0:
        return send_int(len(var))+var.encode()
    else:
        return b''

def recv_int(s):
    return int.from_bytes(s.recv(4), byteorder='little')

def maos_pause(s,pause):
    msg=send_int(21)+send_int(pause)
    s.send(msg)

def maos_get_var(s,name):
    print('receiving ', name)
    msg=send_int(20)+send_int(1) 
    ans=s.send(msg)
    msg=send_str(name)
    ans=s.send(msg)
    ans=readbin.readbin(s) #ans=aolib.readsock(s.fileno())
    print(name, ans.shape, ans.dtype)
    return ans

def maos_get_list(s):
    msg=send_int(20)+send_int(1) 
    ans=s.send(msg)
    msg=send_str("list")
    ans=s.send(msg)
    (ans, header)=readbin.readbin(s, 1) #ans=aolib.readsock(s.fileno())
    ans=header['COMMENT'].rstrip("\n\0").split("\n")[1:]
    return ans

def maos_client(host, port, pid=0):
    s=socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.connect((host, port))
    
    msg=send_int(13)+send_int(pid)
    s.send(msg)
    ans=recv_int(s)

    vars=maos_get_list(s)

    for nf in vars :
        res[nf]=maos_get_var(s,nf)
        
    s.close()
    return res
if __name__ == '__main__':
    if len(sys.argv)>3:
        pid=int(sys.argv[2])
    if len(sys.argv)>2:
        maos_client(sys.argv[1], sys.arg[2], pid)
    else:
        print('Usage: maos_client host port [pid]')
