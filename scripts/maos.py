#!/usr/bin/env python
import re
#parse C structs and its fields, and return variable name and type
def parse_file(srcdir, files):
    lines=list()
    for fn in files:
        fpi=open(srcdir+'/'+fn,'r')
        ln=fpi.read();
        fpi.close();
        while True:
            start=ln.find('#')
            end=ln.find('\n', start)
            if start==-1 or end==-1:
                break
            ln=ln[0:start]+ln[end+1:]

        while True:
            start=ln.find('/*')
            end=ln.find('*/', start)
            if start==-1 or end==-1:
                break
            ln=ln[0:start]+ln[end+2:]
            ln=ln.replace('\n','')

        ln=ln.replace('*','* ') #add space after *
        ln=re.sub(r'[ ]*[*]','*', ln) #remove space before *
        ln=re.sub(r'[ ]+',' ', ln) #remove double space
        lines.append(ln)
    return lines

def parse_structs(srcdir, files):
    lines=parse_file(srcdir, files)
    structs=dict()
    for ln in lines:
        end=0
        while True:
            struct=ln.find('struct', end)
            if struct==-1:
                break
            struct=struct+len('struct ')
            start=ln.find('{', struct)
            if start==-1:
                break;
            name=ln[struct:start];
            end0=ln.find('{', start+1)
            end=ln.find('}', start+1)
            if end0!=-1 and end0<end:
                print(struct, start, end0, end)
                print (ln[struct:end])
                print ("embeded {} are not supported")
                break
            fields=ln[start+1:end].split(';')
            if len(name)==0:
                end0=ln.find(';', end+1);
                name=ln[end+1:end0];
            structs[name]=dict()
            for fi in fields:
                tmp=fi.replace('struct','').replace('const','').lstrip().rstrip().split(' ');
                if len(tmp)==2:
                    structs[name][tmp[1]]=tmp[0]
    return structs

    #parse C functions
def parse_func(srcdir, structs, files):
    funcs=dict()
    lines=parse_file(srcdir, files)
    for ln in lines:
        ln=ln.replace('const','')
        end=0
        while True:
            funname=''
            openp=ln.find('(', end)
            if openp==-1:
                break;
            closep=ln.find(')', openp+1);
            openp2=ln.find('(', openp+1);
            if openp2!=-1 and openp2<closep:
                print(ln[openp:closep])
                print ("embeded () are not supported")
                break
            args=ln[openp+1:closep].split(',')
            defs=ln[end:openp].lstrip().rstrip().split(' ')
            if len(defs)==2:
                funtype=defs[0]
                funname=defs[1]
            else:
                print defs
            funarg=list()
            for arg in args:
                argpair=arg.lstrip().rstrip().split(' ')
                if len(argpair)==2:
                    funarg.append(argpair)
                else:
                    print argpair
            if len(funname)>0:
                funcs[funname]=[funtype, funarg]
            end=closep+1;
            if ln[end]==';' or ln[end]==' ':
                end+=1
    return funcs
            
