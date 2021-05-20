#!/usr/bin/env python3

#Convert functions defined in aolib.h to aolib.c with conversion

from __future__ import print_function
import sys
import maos
from pathlib import Path
if len(sys.argv)>2:
    srcdir=sys.argv[1];
    fnout=sys.argv[2];
else:
    srcdir=str(Path.home())+'/work/programming/aos'
    fnout=srcdir+'/mex/aolib.c'

simu_all=list();

headlist=['maos/parms.h','maos/types.h','lib/accphi.h','lib/cn2est.h','lib/kalman.h',
          'lib/locfft.h','lib/muv.h','lib/servo.h','lib/stfun.h','lib/turbulence.h', 'lib/mkdtf.h','lib/hyst.h']
structs=maos.parse_structs(srcdir, headlist)

funcs=maos.parse_func(srcdir, structs, ['mex/aolib.h'])

fpout=open(fnout,'w')
print('#ifdef __INTEL_COMPILER\n#undef _GNU_SOURCE\n#endif\n#include "interface.h"\n', file=fpout)
def handle_type(argtype):
    if argtype=='dmat*':
        mx2c='mx2d'
        c2mx='any2mx'
        free_c='dfree'
    elif argtype=='cmat*':
        mx2c='mx2c'
        c2mx='any2mx'
        free_c='cfree'
    elif argtype=='lmat*':
        mx2c='mx2l'
        c2mx='any2mx'
        free_c='lfree'
    elif argtype=='dsp*':
        mx2c='mx2dsp'
        c2mx='dsp2mx'
        free_c='dspfree'
    elif argtype=='loc_t*':
        mx2c='mx2loc'
        c2mx='loc2mx'
        free_c='locfree'
    elif argtype=='int' or argtype=='long' or argtype=='double' or argtype=='real' or argtype=='float':
        mx2c='('+argtype+')mxGetScalar'
        c2mx='mxCreateDoubleScalar'
        free_c=''
    elif argtype=='char*':
        mx2c='mx2str'
        c2mx='str2mx'
        free_c='free'
    elif argtype=='dcell*':
        mx2c='mx2dcell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='ccell*':
        mx2c='mx2ccell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='lcell*':
        mx2c='mx2lcell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='loccell*':
        mx2c='mx2loccell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='dspcell*':
        mx2c='mx2dspcell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='cell*':
        mx2c='mx2cell'
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype[-2:]=='**': #output
        mx2c=''
        c2mx='any2mx'
        free_c='cellfree'
    elif argtype=='double*': #output a double (not for vector input)
        mx2c=''
        c2mx='mxCreateDoubleScalar'
        free_c=''
    elif argtype=='rand_t*':
        mx2c='mx2rand'
        c2mx=''
        free_c='free'
    elif argtype=='kalman_t*':
        mx2c='mx2kalman'
        c2mx='kalman2mx'
        free_c='kalman_free'
    elif argtype=='cn2est_t*':
        mx2c='unknown'
        c2mx='cn2est2mx'
        free_c='cn2est_free'
    elif argtype=='dtf_t*':
        mx2c='unknown'
        c2mx='dtf2mx'
        free_c='dtf_free'
    else:
        mx2c='unknown'
        c2mx='unknown'
        free_c=''
    return (mx2c, c2mx,free_c)
funcalls=list()
for funname in funcs: #loop over functions
    funtype=funcs[funname][0]
    funargs=funcs[funname][1]
    funname2=funcs[funname][2]
#    print (funname, funtype,funargs)
    nargs=len(funargs)
    fundef=''
    fundef_free=''
    pointer_output=''
    funout=''
    if funtype=='void':
        pointer_output_count=0
    else:
        pointer_output_count=1
        funout+='out,'

    #handle Input Arguments 
    count=0
    for arg in funargs: #get data from matlab
        argtype=arg[0]
        argname=arg[1]
        mx2c, c2mx, free_c=handle_type(argtype)
        if len(mx2c)>0: #input argument
            fundef+='    '+argtype+' '+argname+'='+mx2c+'(prhs['+str(count)+']);\n' #input from matlab
            count=count+1
        else: #output argument
            fundef+='    '+argtype[0:-1]+' '+argname+'=0;\n' #output
            arg[1]='&'+argname
            pointer_output+='    plhs['+str(pointer_output_count)+']='+c2mx+'('+argname+');\n'
            funout+=argname+','
            pointer_output_count+=1
            nargs-=1;
        if len(free_c)>0:
            fundef_free+='    '+free_c+'('+argname+');\n'

    #Call the C function
    if funtype=='void':
        fundef+='    '+funname2+"("
    else:
        fundef+='    '+funtype+' '+funname+'_out='+funname2+"("
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
    fundef0+='    if(nrhs!='+str(nargs)+') {\n'
    funcall="    printf(\""
    if len(funout)>0:
        funcall+="["+funout[0:-1]+"]="

    funcall+="aolib('"+funname+"',"
    for arg in funargs:
        argname=arg[1]
        if argname[0]!='&':
            funcall+=argname+","
    funcall=funcall[0:-1]+")\\n\");\n"  
    funcalls.append(funcall);
    fundef0+="    "+funcall;
    fundef0+='        mexErrMsgTxt(\"Expect '+str(nargs)+' arguments.");\n    }\n'
    fundef=fundef0+fundef
    print(fundef, file=fpout)

print("void print_usage(){\n    printf(\"Usage:\\n\");", file=fpout)
funcalls.sort()
print('\n'.join(funcalls), file=fpout)
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
