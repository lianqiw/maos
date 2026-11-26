#!/usr/bin/env python3
"""Generate mex/maos2mex.h, which exports the maos simulation state to Matlab.

Usage:
    maos2mex.py [srcdir output]

    srcdir: MAOS source tree (the directory containing maos/, lib/, mex/).
    output: file to write.
    Without arguments, srcdir is the parent directory of this script and the
    output is maos2mex.h in the current directory.

    The build runs it from the mex build directory (mex/Makefile.am) as
        scripts/maos2mex.py $(SRC_DIR) maos2mex.h

Inputs:
    Struct definitions parsed (with maos_parse.py) from srcdir/lib/*.h,
    sys/scheduler_client.h, mex/aolib.h, maos/parms.h and maos/types.h.
    Starting from sim_t, struct members that are themselves parsed structs
    are expanded recursively.

Output:
    C code included by the maos mex (mex/maos.c) containing, for each
    reachable struct type, "static mxArray *get_<type>(...)" which converts
    the struct to a Matlab struct (members of supported types only; others
    are skipped), get_data_help(), and get_data(simu, key) which returns the
    struct selected by key (e.g. "parms", "recon"), as used by
    maos('get', key) in Matlab.
"""
from __future__ import print_function
import sys,os.path
import glob
import maos_parse
from pathlib import Path
if len(sys.argv)==3:
    srcdir=sys.argv[1]
    fnout=sys.argv[2]
elif len(sys.argv)==1:
    srcdir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fnout='maos2mex.h' #output to build directory
else:
    print('Usage: '+sys.argv[0]+' [srcdir output]', file=sys.stderr)
    sys.exit(0 if len(sys.argv)==2 and sys.argv[1] in ('-h','--help') else 1)
if not os.path.isdir(srcdir+'/maos/'):
    sys.exit('maos2mex.py: '+srcdir+' is not the maos source dir (no maos/ subdirectory)')

simu_all=list();

headerlist=glob.glob(srcdir+'/lib/*.h')
headerlist.append(srcdir+'/sys/scheduler_client.h')
headerlist.append(srcdir+'/mex/aolib.h')
headerlist.append(srcdir+'/maos/parms.h')
headerlist.append(srcdir+'/maos/types.h')
structs=maos_parse.parse_structs('', headerlist)
#funcs=maos_parse.parse_func('',headerlist)

if 'sim_t' not in structs:
    sys.exit('maos2mex.py: unable to parse struct sim_t from '+srcdir+'/maos/types.h')
simu=structs['sim_t']
unknowns=dict() #count of unsupported member types (for debugging)
expanding=set() #id() of structs being expanded, to stop recursive types
def expand_struct(struct):
    """Expand, in place, struct members whose type is a parsed struct.

    Parameters:
        struct: dict {member: type} from maos_parse.parse_structs(). Each
            member whose type (with one '*' removed) is a parsed struct is
            replaced by [type, member dict], recursively. Members already
            expanded (lists) are left alone. A member of a struct type that
            is currently being expanded (a recursive type) is left as a type
            string, so it is skipped as an unknown pointer by var2mx().
    """
    expanding.add(id(struct))
    for key in struct:
        val0=struct[key]
    
        if isinstance(val0, str) and val0[-1]=='*':
            val=val0[0:-1]
            ispointer=1
        else:
            val=val0
            ispointer=0
        if isinstance(val, str) and structs.get(val):
            if id(structs[val]) in expanding: #recursive type. Do not expand.
                continue
            if isinstance(structs[val], dict): #other types
                expand_struct(structs[val])
            struct[key]=[val0, structs[val]]
    expanding.discard(id(struct))

expand_struct(simu)
#Convert cname of type ctype in C to mexname in Matlab
def var2mx(mexname, cname, ctype):
    """Return a C statement converting a non-struct member to an mxArray.

    Parameters:
        mexname: C lvalue receiving the mxArray (e.g. "\ttmp2").
        cname: C expression of the member (e.g. "simu->dmreal").
        ctype: its type as parsed (e.g. "dcell*", "int").

    Returns:
        "mexname=xxx2mx(cname);" for numeric arrays/cells/loc/map/sparse
        (any2mx), strings (str2mx), and mxCreateDoubleScalar() for scalars.
        For other types a comment starting with "//" is returned (the caller
        drops it) and the type is counted in unknowns.
    """
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
            if unknowns.get(ctype0):
                unknowns[ctype0]+=1
            else:
                unknowns[ctype0]=1
            return "//unknown pointer "+cname+":"+ctype0+';'
    else:
        if ctype=='char' or ctype=='int' or ctype=='long' or ctype=='double' or ctype=='float' or ctype=='real':
            return mexname+'=mxCreateDoubleScalar('+cname+');'
        else:
            if unknowns.get(ctype0):
                unknowns[ctype0]+=1
            else:
                unknowns[ctype0]=1
            return "//unknown type    "+cname+":"+ctype0+';'

fundefs=dict()    #get_<type> name -> function definition
funheaders=dict() #get_<type> name -> prototype
funcalls=dict()   #get_<type> name -> call used by get_data() for top level members
def struct2mx(vartype, nvar, funprefix, parentpointer, fullpointer, varname, struct):
    """Generate get_<type>() converting a struct to a Matlab struct.

    The definition and prototype are stored in fundefs and funheaders (keyed
    by the type name without "_t" and "*", lower case; a type reached again
    overwrites the same entry). Members that are expanded structs are
    handled recursively; the member "save" is skipped. Members named powfs,
    wfs and wfsr are arrays of structs whose length is (parms->)npowfs,
    nwfs or nwfsr; all others are single structs.

    Parameters:
        vartype: C type of the struct, with or without '*'. Without '*' the
            struct is embedded and its address is passed.
        nvar: ignored; recomputed from varname as described above.
        funprefix: "" for the top level call, otherwise the parent's type
            name plus "_"; used to decide which members get_data() exports.
        parentpointer: C expression prefix to reach this member from the
            parent's get_ function (e.g. "simu->").
        fullpointer: C expression prefix to reach this member from simu.
        varname: member name (e.g. "parms").
        struct: member dict from expand_struct().

    Returns:
        The C call of the generated function (e.g. "get_parms(simu->parms);"),
        or '' if the struct has no member that could be converted.
    """
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
    nvalid=0
    for key in struct: 
        valtype=struct[key]
        ans=""
        if isinstance(valtype, list):
            if key=='save':
                continue
            tmp=struct2mx(valtype[0], nvar, funname+"_", childpointer, fullpointer+childpointer, key, valtype[1])
            if len(tmp)>0 and tmp[0:2]!='//':
                ans="\ttmp2="+tmp
            else:
                ans="//unable to parse struct"+valtype[0]
        else:
            ans=var2mx("\ttmp2", childpointer+key, valtype)
        if len(ans)>0 and ans[0:2]!='//':
            funbody.append(ans+"i=mxAddField(tmp, \""+key+"\");mxSetFieldByNumber(tmp, "+stind+", i, tmp2);")
            nvalid+=1
    funbody.sort()
    fundef+="\n".join(funbody)
    if nvar!="1":
        fundef+="\n}\n"
    fundef+="\n\treturn tmp;\n}\n"
    if nvalid==0:
        return ''
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
#for x in unknowns:
    #print("  Unknown ({:2d}) times type: {}".format(unknowns[x], x))

#Sort key: names starting with a number sort numerically by their first
#space separated word, before all others; then alphabetically.
def sort_key(item):
    return (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)
with open(fnout,'w') as fc:
    keys=funheaders.keys()
    keys=sorted(keys, key=sort_key)
    for key in keys:
        print(funheaders[key], file=fc)
    keys=fundefs.keys()
    keys=sorted(keys, key=sort_key)
    for key in keys:
        print(fundefs[key], file=fc)
    print("static void get_data_help(){\n", end="", file=fc)
    for key in funcalls:
        print("\tinfo(\""+key+"=maos('get','"+key+"')\\n\");", file=fc)
    print("}", file=fc)
    print("static mxArray *get_data(sim_t *simu, const char *key){\n\t", end="", file=fc)
    for key in funcalls:
        #print(key)
        print("if(!strcmp(key, \""+key.replace("simu_","")+"\")) return "+funcalls[key]+"\n\telse ", end="", file=fc)
    print("{get_data_help();return mxCreateDoubleMatrix(0,0,mxREAL);}\n}", file=fc)
