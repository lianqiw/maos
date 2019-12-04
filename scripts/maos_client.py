#!env python3
import aolib
import readbin
#get variables from running MAOS

def send_int(var):
    return var.to_bytes(4, byteorder='little')
def send_str(var):
    if len(var)>0:
        return send_int(len(var))+var.encode()
    else:
        return b''

import socket
s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(('maxwell', 10000))
pid=0
msg=send_int(13)+send_int(pid)
s.send(msg)
ans=s.recv(4)
print("recv",int.from_bytes(ans, byteorder='little'))
msg=send_int(20)+send_int(1) 
ans=s.send(msg)
print("send",ans)
msg=send_str("dmreal")
ans=s.send(msg)
print("send",ans)
#ans=s.recv(4)
#print("recv",int.from_bytes(ans, byteorder='little'))
#dm=readbin.readbin(s)
dm=aolib.readsock(s.fileno())
dm2=dm*2;
msg=send_int(20)+send_int(2) 
ans=s.send(msg)
msg=send_str("dmreal")
ans=s.send(msg)
aolib.writesock(dm2, s.fileno())


s.close()
