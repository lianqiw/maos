#!/usr/bin/env python
import re
import sys
fp=list();
simu_all=list();
fp.append(open('/Users/lianqiw/work/programming/aos/maos/parms.h','r'));
fp.append(open('/Users/lianqiw/work/programming/aos/maos/types.h','r'));
structs=dict()
for fpi in fp:
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
    end=0
    while True:
        struct=ln.find('struct', end)
        if struct==-1:
            break
        struct=struct+len('struct ')
        start=ln.find('{', struct)
        name=ln[struct:start];
        end=ln.find('}', start)
        structs[name]=dict()
        fields=ln[start+1:end].split(';')
        for fi in fields:
            tmp=fi.replace('struct','').replace('const','').lstrip().rstrip().split(' ');
            if len(tmp)==2:
                structs[name][tmp[1]]=tmp[0]
simu=dict()
simu_type=dict() 
    
for key in structs['SIM_T']:
    val=structs['SIM_T'][key]
    if val[-1]=='*':
        val=val[0:-1];
        ispointer=1
    else:
        ispointer=0
    if structs.has_key(val): #other types
        simu[key]=structs[val]
        if ispointer==1:
           simu_type[key]='->'
        else:
           simu_type[key]='.'
    else:
        if ispointer==1:
            simu[key]=val+'*'
            simu_type[key]='*'
        else:
            simu_type[key]=''
            simu[key]=val
if len(sys.argv)>1:
    fnout=sys.argv[1]
else:
    fnout='typeconv.h'
fc=open(fnout,'w')
def var2mx(mexname, cname, ctype):
    ans=1
    if ctype[-1]=='*':
        ccast=''
        ctype=ctype[0:-1]
        if ctype[-2:]=='_t':
            ctype=ctype[0:-2];
        if ctype=='map':
            ctype='dmat'
            ccast='(const dmat*)'
        if ctype=='dmat' or ctype=='cmat':
            fun_c=ctype[0]
        elif ctype=='loc' or ctype=='dcell' or ctype=='ccell':
            fun_c=ctype
        else:
            fun_c=''
        if len(fun_c)>0:
            print>>fc, mexname+"="+fun_c+'2mx('+ccast+cname+');'
        else:
            print>>fc, "//unknown1 type "+cname+":"+ctype
            ans=0
    else:
        if ctype=='int' or ctype=='double':
            print>>fc, mexname+'=mxCreateDoubleScalar('+cname+');'
        else:
            print>>fc, "//unknown2 type "+cname+":"+ctype
            ans=0
    return ans
print>>fc, "mxArray *get_simu(SIM_T *simu){"
print>>fc, "if(!simu) return mxCreateDoubleMatrix(0,0,mxREAL);"
print>>fc, "int i;\nint i_n;"
print>>fc, "mxArray *tmp=0, *tmp2=0;"
print>>fc, "int npowfs=simu->parms->npowfs;"
print>>fc, "mxArray *out=mxCreateStructMatrix(1,1,0,0);"
for key in simu:
    val=simu[key]
    count=0
    if type(val)==type(dict()):
        if key=='save':
            continue
        if simu_type[key]=='->':
            if key=='powfs':
                n_name="npowfs"
                n=2
            else:
                n_name="1"
                n=1
        else:
            n=-1
        print>>fc, "tmp=mxCreateStructMatrix("+str(n)+",1,0,0);"
        for key2 in val:
            val2=val[key2]
            if n>0:
                if n>1:
                    print>>fc, "for(i_n=0; i_n<"+n_name+"; i_n++){"
                    key2full='simu->'+key+'[i_n].'+key2
                else:
                    key2full='simu->'+key+'->'+key2
            else:
                key2full='simu->'+key+'.'+key2

            print key2full
            if var2mx("tmp2", key2full, val2):
                print>>fc, "i=mxAddField(tmp, \""+key2+"\");"
                if n>1:
                    print>>fc, "mxSetFieldByNumber(tmp, i_n, i, tmp2);"
                else:
                    print>>fc, "mxSetFieldByNumber(tmp, 0, i, tmp2);"
                count=count+1
            if n>1:
                print >>fc, "}"
    else:
        keyfull='simu->'+key
        print keyfull
        if var2mx("tmp", keyfull, val):
            count=count+1
    if count>0:
        print>>fc, "i=mxAddField(out, \""+key+"\");"
        print>>fc, "mxSetFieldByNumber(out, 0, i, tmp);"

print>>fc, "return out;\n}\n"
fc.close()
