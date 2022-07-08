#!/usr/bin/env python
import re

def parse_file(srcdir, files):
    lines=list()
    for fn in files:
        fpi=open(srcdir+'/'+fn,'r')
        ln=fpi.read();
        fpi.close();
    
        while True: #remove all comments start with # or //
            start=ln.find('#include')
            if start==-1:
                start=ln.find('//')
            if start==-1:
                start=ln.find('#define')
            if start==-1:
                break
            end=ln.find('\n', start)
            if end==-1:
                break
            ln=ln[0:start]+ln[end+1:]

        while True: #Remove all C style comments.
            start=ln.find('/*')
            end=ln.find('*/', start)
            if start==-1 or end==-1:
                break
            ln=ln[0:start]+ln[end+2:]

        ln=ln.replace('\n',' ')
        ln=ln.replace(';','; ')
        ln=re.sub(r'format,[ ]*...','format',ln)
        ln=re.sub(r'[)][ ]*[(]',')(',ln) #replace space between ) and (
        ln=ln.replace('*','* ') #add space after *
        ln=re.sub(r'[ ]*[*]','*', ln) #remove space before *
        ln=re.sub(r'[ ]+',' ', ln) #remove double space
        lines.append(ln)
    return lines
#parse C structs and its fields, and return variable name and type
def parse_structs(srcdir, files):
    lines=parse_file(srcdir, files)
    structs=dict()
    for ln in lines:
        end=0
        while True:
            struct=ln.find('struct', end)
            if struct==-1:
                break
            struct=struct+len('struct')
            start=ln.find('{', struct)
            if start==-1:
                break;
            name=ln[struct+1:start];
            end0=ln.find('{', start+1)
            end=ln.find('}', start+1)

            if end0!=-1 and end0<end:
                #print(struct, start, end0, end)
                #print (ln[struct:end])
                #print ("embeded {} are not supported")
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
        if ln[0:6]=='extern':
            pass
        ln=ln.replace('const','')
        end=0
        while True:
            funname=''
            openp=ln.find('(', end)
            if openp==-1:
                break;
            closep=ln.find(')', openp+1);
            openp2=ln.find('(', openp+1);
            if openp2!=-1 and openp2<closep: #nested (). Skip
                closep2=ln.find(')', closep+1)
                end=closep2+1
                continue
            elif openp2==closep+1: #()()
                closep2=ln.find(')', openp2+2)
                end=closep2+1
                continue

            args=ln[openp+1:closep].split(',')
            defs=ln[end:openp].lstrip().rstrip().split(' ')
            msg=''
            if len(defs)>=2:
                funtype=defs[-2]
                funname=defs[-1]
                name2=funname.split('=')
                funname=name2[0] #Mex name
                funname2=name2[len(name2)-1] #C name
            else:
                msg='ill function definition('+ln[end:openp]+')'
                funname=''
            if funname=='readbin':
                print(args)
            funarg=list()
            for arg in args:
                argpair=arg.lstrip().rstrip().split(' ')
                if len(argpair)==1:
                    if argpair[0]=='...':
                        pass #skip this parameter
                    else:
                        try:
                            junk=float(argpair[0])
                            funarg.append(['number', argpair[0]])
                        except:
                            msg='ill argument('+arg+')'
                            funname=''
                elif len(argpair)==2 :
                    funarg.append(argpair)
                else:#failed to parse function
                    msg='ill argument('+arg+')'
                    raise(Exception(msg))
                    funname=''
            if funname=='readbin':
                print(args)
            if len(funname)>0:
                funcs[funname]=[funtype, funarg, funname2]
            else:
                print("Skipped:", ln[end:closep], msg,"\n")
            end=closep+1;
            if ln[end]==';' or ln[end]==' ':
                end+=1
    return funcs
            
