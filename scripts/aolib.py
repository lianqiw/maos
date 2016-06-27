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
def handle_type(argtype):
    if argtype=='dmat*':
        mx2c='mx2d'
        c2mx='any2mx'
        free_c='dfree'
    elif argtype=='dsp*':
        mx2c='mx2dsp'
        c2mx='dsp2mx'
        free_c='dspfree'
    elif argtype=='loc_t*':
        mx2c='mx2loc'
        c2mx='loc2mx'
        free_c='locfree'
    elif argtype=='int' or argtype=='long' or argtype=='double':
        mx2c='('+argtype+')mxGetScalar'
        c2mx='mxCreateDoubleScalar'
        free_c=''
    elif argtype=='char*':
        mx2c='mx2str'
        c2mx='str2mx'
        free_c='free'
    elif argtype=='void*' or argtype=='cell*' or argtype=='dcell*' or argtype=='ccell*':
        mx2c='mx2any'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='cell**':
        mx2c=''
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='rand_t*':
        mx2c='mx2rand'
        c2mx=''
        free_c='free'
    else:
        mx2c='unknown'
        c2mx='unknown'
        free_c=''
    return (mx2c, c2mx,free_c)

for funname in funcs: #loop over functions
    funtype=funcs[funname][0]
    funargs=funcs[funname][1]
    print (funname, funtype,funargs)
    nargs=len(funargs)
    fundef=''
    fundef_free=''
    pointer_output=''
    if funtype=='void':
        pointer_output_count=0
    else:
        pointer_output_count=1
        
    count=0
    for arg in funargs: #get data from matlab
        argtype=arg[0]
        argname=arg[1]
        mx2c, c2mx, free_c=handle_type(argtype)
        if len(mx2c)>0:
            fundef+='    '+argtype+' '+argname+'='+mx2c+'(prhs['+str(count)+']);\n' #input frm matlab
            count=count+1
        else:
            fundef+='    '+argtype[0:-1]+' '+argname+'=0;\n' #output
            arg[1]='&'+argname
            pointer_output+='    plhs['+str(pointer_output_count)+']='+c2mx+'('+argname+');\n'
            pointer_output_count+=1
            nargs-=1;
        if len(free_c)>0:
            fundef_free+='    '+free_c+'('+argname+');\n'
    if funtype=='void':
        fundef+='    '+funname+"("
    else:
        fundef+='    '+funtype+' '+funname+'_out='+funname+"("
    for arg in funargs:
        argname=arg[1]
        fundef+=argname+","
    fundef=fundef[0:-1]+');\n'
    mx2c, c2mx, free_c=handle_type(funtype)
    if funtype !='void':
        fundef+='    plhs[0]='+c2mx+'('+funname+'_out);\n'
    else:
        fundef+='    (void)plhs;\n'
    if len(pointer_output)>0:
        fundef+=pointer_output
    if len(free_c)>0:
        fundef+='    '+free_c+'('+funname+'_out);\n'
    fundef+=fundef_free
    fundef+='}'
    fundef0='void '+funname+'_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){\n'
    fundef0+='    (void)nlhs;\n';
    fundef0+='    if(nrhs!='+str(nargs)+') mexErrMsgTxt(\"Expect '+str(nargs)+' arguments\\n");\n'
    fundef=fundef0+fundef
    print(fundef, file=fpout)

print("void print_usage(){\n    printf(\"Usage:\\n\");", file=fpout)
for funname in funcs:
    funtype=funcs[funname][0]
    funargs=funcs[funname][1]
    if funtype=='void':
        funcall="    printf(\"    aolib('"+funname+"',"        
    else:
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
