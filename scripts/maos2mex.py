#!/usr/bin/env python3
from __future__ import print_function
import sys
import maos
if len(sys.argv)>2:
    srcdir=sys.argv[1];
    fnout=sys.argv[2];
else:
    srcdir='/Users/lianqiw/work/programming/aos'
    fnout='../mex/maos2mex.h'

simu_all=list();

headlist=['maos/parms.h','maos/types.h','lib/accphi.h','lib/cn2est.h','lib/kalman.h',
          'lib/locfft.h','lib/muv.h','lib/servo.h','lib/stfun.h','lib/turbulence.h',
          'lib/mkdtf.h', 'math/chol.h','sys/scheduler.h']

#Obtain the definition of all structs
dictall=maos.parse_structs(srcdir, headlist)

simu=dict()
simu=dictall['sim_t']

def expand_struct(struct):
    for key in struct:
        val0=struct[key]
    
        if type(val0)==type('') and val0[-1]=='*':
            val=val0[0:-1]
            ispointer=1
        else:
            val=val0
            ispointer=0
        if type(val)==type('') and dictall.get(val):
            if type(dictall[val])==type(dict()): #other types
                expand_struct(dictall[val])
            struct[key]=[val0, dictall[val]]

expand_struct(simu)
#Convert cname of type ctype in C to mexname in Matlab
def var2mx(mexname, cname, ctype):
    out=""
    ans=1
    ctype0=ctype
    if ctype[-1]=='*':
        ccast=''
        ctype=ctype[0:-1]
        if ctype[-2:]=='_t':
            ctype=ctype[0:-2];
        if ctype=='map' or ctype[1:]=='mat' or ctype[-4:]=='cell' or ctype=='loc' or ctype=='pts' or ctype[1:]=='sp':
            fun_c='any'
        elif ctype=='char':
            fun_c='str'
        else:
            fun_c=''
        if len(fun_c)>0:
            return mexname+"="+fun_c+'2mx('+ccast+cname+');'
        else:
            print("//unknown pointer "+cname+":"+ctype0)
            return "//unknown pointer "+cname+":"+ctype0+';'
    else:
        if ctype=='char' or ctype=='int' or ctype=='long' or ctype=='double' or ctype=='float' or ctype=='real':
            return mexname+'=mxCreateDoubleScalar('+cname+');'
        else:
            print("//unknown type    "+cname+":"+ctype0)
            return "//unknown type    "+cname+":"+ctype0+';'

fundefs=dict()
funheaders=dict()
funcalls=dict()
def struct2mx(vartype, nvar, funprefix, parentpointer, fullpointer, varname, struct):
    if vartype[-1]=='*':
        ispointer=1
    else:
        ispointer=0
        vartype=vartype+"*"

    if varname=="powfs":
        if parentpointer=="simu->":
            nvar="parms->npowfs"
        else:
            nvar="npowfs"
    elif varname=="wfs":
        if parentpointer=="simu->":
            nvar="parms->nwfs"
        else:
            nvar="nwfs"
    elif varname=="wfsr":
        if parentpointer=="simu->":
            nvar="parms->nwfsr"
        else:
            nvar="nwfsr"
    else:
        nvar="1"
    if nvar=="1":
        childpointer=varname+'->'
    else:
        childpointer=varname+'[i_n].'
#    if funprefix=="simu_":
#        funname=varname
#    else:
#        funname=funprefix+varname;
    funname=vartype.replace('_t','').replace('*','').lower();
    if nvar!="1":
        funheader="static mxArray *get_"+funname+"(const "+vartype+" "+varname+", int nvar);";
    else:
        funheader="static mxArray *get_"+funname+"(const "+vartype+" "+varname+");";
    fundef=funheader[0:-1]+"{\n"
    fundef+="\tif(!"+varname+") return mxCreateDoubleMatrix(0,0,mxREAL);\n"
    if nvar!="1":
        fundef+="\tmxArray* tmp=mxCreateStructMatrix(nvar,1,0,0);\n\tint i;\n\tmxArray *tmp2;\n"
    else:
        fundef+="\tmxArray* tmp=mxCreateStructMatrix(1,1,0,0);\n\tint i;\n\tmxArray *tmp2;\n"
    if nvar!="1":
        fundef+="for(int i_n=0; i_n<nvar; i_n++){\n"
        stind='i_n'
    else:
        stind='0'
    funbody=list()
    for key in struct: 
        valtype=struct[key]
        ans=""
        if type(valtype)==type(list()):
            if key=='save':
                continue
            ans+="\ttmp2="+struct2mx(valtype[0], nvar, funname+"_", childpointer, fullpointer+childpointer, key, valtype[1])
        else:
            ans+=var2mx("\ttmp2", childpointer+key, valtype)
        if len(ans)>0:
            funbody.append(ans+"i=mxAddField(tmp, \""+key+"\");mxSetFieldByNumber(tmp, "+stind+", i, tmp2);")
    funbody.sort()
    fundef+="\n".join(funbody)
    if nvar!="1":
        fundef+="\n}\n"
    fundef+="\n\treturn tmp;\n}\n"
    fundefs[funname]=fundef
    funheaders[funname]=funheader
    #print(funname, funprefix)
    #if funname.replace('simu_','').find('_')==-1:
    if (funprefix=='' or funprefix=="sim_") and varname.find('_')==-1:
        if ispointer==0:
            funcalls[funname]="get_"+funname+"(&("+fullpointer+varname+"));"
        elif nvar!="1":
            funcalls[funname]="get_"+funname+"("+fullpointer+varname+", "+fullpointer+nvar+");"
        else:
            funcalls[funname]="get_"+funname+"("+fullpointer+varname+");"

    if ispointer==0:
        return "get_"+funname+"(&("+parentpointer+varname+"));"
    elif nvar!="1":
        return "get_"+funname+"("+parentpointer+varname+", "+parentpointer+nvar+");"
    else:
        return "get_"+funname+"("+parentpointer+varname+");"

struct2mx("sim_t*", "1", "","","", "simu", simu)


fc=open(fnout,'w')
keys=funheaders.keys()
keys=sorted(keys, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))
for key in keys:
    print(funheaders[key], file=fc)
keys=fundefs.keys()
keys=sorted(keys, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item))
for key in keys:
    print(fundefs[key], file=fc)
print("static void get_data_help(){\n", end="", file=fc)
for key in funcalls:
    print("\tinfo(\""+key+"=maos('get','"+key+"')\\n\");", file=fc)
print("}", file=fc)
print("static mxArray *get_data(sim_t *simu, const char *key){\n\t", end="", file=fc)
for key in funcalls:
    print(key)
    print("if(!strcmp(key, \""+key.replace("simu_","")+"\")) return "+funcalls[key]+"\n\telse ", end="", file=fc)
print("{get_data_help();return mxCreateDoubleMatrix(0,0,mxREAL);}\n}", file=fc)
fc.close()
