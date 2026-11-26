#!/usr/bin/env python3
"""Generate mex/aolib.c, the Matlab MEX gateway of the aolib functions.

Usage:
    lib2mex.py [srcdir output]

    srcdir: MAOS source tree (the directory containing maos/, lib/, mex/).
    output: file to write.
    Without arguments, srcdir is the parent directory of this script and the
    output is aolib.c in the current directory.

    The build runs it from the mex build directory (mex/Makefile.am) as
        scripts/lib2mex.py $(SRC_DIR) aolib.c

Inputs:
    The functions to export are parsed (with maos_parse.py) from
    srcdir/mex/aolib.h, where "type mexname=cname(args);" exports C function
    cname as mexname. Struct definitions are also parsed from lib/*.h,
    sys/scheduler_client.h, mex/aolib.h, maos/parms.h and maos/types.h, but
    are currently not used.

Output:
    C source that includes interface.h and defines, for each function,
    <name>_mex(nlhs, plhs, nrhs, prhs) which converts the Matlab inputs with
    mx2* helpers, calls the C function, converts the result and any
    pointer-to-pointer (type**) or double* outputs to plhs[] (C return value
    first), and frees temporaries. It also defines print_usage() and
    mexFunction(), which dispatches on the first (string) argument, so the
    Matlab call is aolib('name', args...). The file is written only after it
    has been generated completely.

An argument type without a mapping in handle_type() produces C code that
calls unknown() and does not compile; a warning is printed to stderr.
"""

from __future__ import print_function
import sys,os.path
import io
import glob
import maos_parse
from pathlib import Path
if len(sys.argv)==3:
    srcdir=sys.argv[1]
    fnout=sys.argv[2]
elif len(sys.argv)==1:
    srcdir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fnout='aolib.c'
else:
    print('Usage: '+sys.argv[0]+' [srcdir output]', file=sys.stderr)
    sys.exit(0 if len(sys.argv)==2 and sys.argv[1] in ('-h','--help') else 1)
if not os.path.isdir(srcdir+'/maos/'):
    sys.exit('lib2mex.py: '+srcdir+' is not the maos source dir (no maos/ subdirectory)')
if not os.path.isfile(srcdir+'/mex/aolib.h'):
    sys.exit('lib2mex.py: '+srcdir+'/mex/aolib.h not found')
simu_all=list()
headerlist=glob.glob(srcdir+'/lib/*.h')
headerlist.append(srcdir+'/sys/scheduler_client.h')
headerlist.append(srcdir+'/mex/aolib.h')
headerlist.append(srcdir+'/maos/parms.h')
headerlist.append(srcdir+'/maos/types.h')
structs=maos_parse.parse_structs('', headerlist)
funcs=maos_parse.parse_func(srcdir,['mex/aolib.h'])

fpout=io.StringIO() #written to fnout at the end, so a failure leaves no partial file
print('#ifdef __INTEL_COMPILER\n#undef _GNU_SOURCE\n#endif\n#include "interface.h"\n', file=fpout)
def handle_type(argtype):
    """Return the conversion functions for a C type.

    Parameters:
        argtype: C type as returned by maos_parse.parse_func (e.g. "dmat*").

    Returns:
        (mx2c, c2mx, free_c): names of the function converting an mxArray to
        the C type (empty for an output argument, i.e. type** or double*),
        the function converting the C value to an mxArray, and the function
        freeing the C value (empty if nothing to free). 'unknown' marks a
        missing conversion.
    """
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
funcalls=list() #usage lines for print_usage()
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
        if argtype=='number': #literal number: passed to C, not taken from Matlab
            nargs-=1
            continue
        mx2c, c2mx, free_c=handle_type(argtype)
        if mx2c=='unknown':
            print('lib2mex.py: warning: no Matlab conversion for argument '+argtype+' '+argname+' of '+funname, file=sys.stderr)
        if len(mx2c)>0: #input argument
            fundef+='    '+argtype+' '+argname+'='+mx2c+'(prhs['+str(count)+']);\n' #input from matlab
            count=count+1
        else: #output argument
            fundef+='    '+argtype[0:-1]+' '+argname+'=0;\n' #output
            arg[1]='&'+argname #pass by address in the C call below
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
    if funtype!='void' and c2mx=='unknown':
        print('lib2mex.py: warning: no Matlab conversion for return type '+funtype+' of '+funname, file=sys.stderr)
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
    for arg in funargs: #Matlab inputs only: skip outputs (&name) and literal numbers
        argname=arg[1]
        if argname[0]!='&' and arg[0]!='number':
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
with open(fnout,'w') as fp:
    fp.write(fpout.getvalue())
fpout.close()
