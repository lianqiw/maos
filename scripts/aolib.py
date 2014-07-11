#!/usr/bin/env python
from __future__ import print_function
import sys
import maos
if len(sys.argv)>2:
    srcdir=sys.argv[1];
    fnout=sys.argv[2];
else:
    srcdir='/Users/lianqiw/work/programming/aos'
    fnout='../mex/aolib.c'

simu_all=list();

headlist=['maos/parms.h','maos/types.h','lib/accphi.h','lib/cn2est.h','lib/kalman.h',
          'lib/locfft.h','lib/muv.h','lib/servo.h','lib/stfun.h','lib/turbulence.h']
structs=maos.parse_structs(srcdir, headlist)

funcs=maos.parse_func(srcdir, structs, ['mex/aolib.h'])

fpout=open(fnout,'w')
print('#ifdef __INTEL_COMPILER\n#undef _GNU_SOURCE\n#endif\n#include "interface.h"\n', file=fpout)

for funname in funcs:
    funtype=funcs[funname][0]
    funargs=funcs[funname][1]
    print (funname, funtype,funargs)
    fundef='void '+funname+'_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){\n'
    fundef+='    if(nrhs!='+str(len(funargs))+') mexErrMsgTxt(\"Expect '+str(len(funargs))+' arguments\\n");\n'
    count=0
    for arg in funargs:
        argtype=arg[0]
        argname=arg[1]
        if argtype=='dmat*':
            fundef+='    '+argtype+' '+argname+'=mx2d(prhs['+str(count)+']);\n'
        count=count+1
    fundef+='    '+funtype+' '+funname+'_out='+funname+"("
    for arg in funargs:
        argname=arg[1]
        fundef+=argname+","
    fundef=fundef[0:-1]+');\n'
    fundef+='    plhs[0]=any2mx('+funname+'_out);\n'
    for arg in funargs:
        argtype=arg[0]
        argname=arg[1]
        if argtype=='dmat*':
            fundef+='    dfree('+argname+');\n'
    fundef+='}'
    print(fundef, file=fpout)

print("void print_usage(){\n    printf(\"Usage:\\n\");", file=fpout)
for funname in funcs:
    funtype=funcs[funname][0]
    funargs=funcs[funname][1]
    funcall="    printf(\"out=aolib('"+funname+"',"
    for arg in funargs:
        argtype=arg[0]
        argname=arg[1]
        funcall+=argname+","
    funcall=funcall[0:-1]+")\\n\");"
    print(funcall, file=fpout)
print("}", file=fpout)
print('''void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs<1){
	print_usage();
        return;
    }
    char *cmd=mxArrayToString(prhs[0]);\n    ''',end="", file=fpout)
for funname in funcs:
    print("if(!strcmp(cmd, \""+funname+"\")) "+funname+"_mex(nlhs, plhs, nrhs-1, prhs+1);\n    else ", end="",file=fpout)
print("print_usage();", file=fpout)
print("}\n", file=fpout)
fpout.close()
