#!/usr/bin/env python3
"""Minimal C header parser shared by the MAOS code generators.

This module is imported (not run) by lib2py.py, lib2mex.py and maos2mex.py,
which use it to extract struct and function declarations from MAOS headers
(lib/*.h, sys/scheduler_client.h, mex/aolib.h, maos/parms.h, maos/types.h).

It is not a real C parser. Each header is flattened to a single line with
comments, #include and #define lines removed, and declarations are then found
with plain string searches. It relies on the formatting conventions of the
MAOS headers:

- struct bodies are not nested. Parsing of a file stops at the first struct
  whose body contains another '{' (e.g. an anonymous union), so structs that
  follow it in the same file are not seen.
- Every function argument has the form "type name" (after removing const and
  joining '*' to the type), or is a literal number or "...".
- mex/aolib.h may use "mexname=cname(args);" to export the C function cname
  under a different name.

Functions:
    parse_file(srcdir, files): read and normalize header text.
    parse_structs(srcdir, files): return {struct name: {field: type}}.
    parse_func(srcdir, files): return {name: [return type, args, C name]}.
"""
import re

def parse_file(srcdir, files):
    """Read C headers and return their text normalized for simple parsing.

    Parameters:
        srcdir: directory prepended (with '/') to each entry in files, or
            None/'' to use the entries as given.
        files: list of header file names.

    Returns:
        A list with one string per file, in the same order. Each string is
        the whole file on a single line with #include/#define lines, // and
        C comments removed, whitespace collapsed to single spaces, and '*'
        attached to the preceding type (e.g. "const dmat *a" -> "const dmat* a").
        "format, ..." is reduced to "format" to drop variadic arguments.

    Raises:
        OSError: if a file cannot be read.
    """
    lines=list()
    for fn in files:
        if srcdir is not None and len(srcdir)>0:
            fn=srcdir+'/'+fn
        with open(fn,'r') as fpi:
            ln=fpi.read()
    
        #Remove lines (from the marker to end of line) starting with #include,
        #// or #define. Only the first marker type still present is searched
        #in each pass, so this repeats until none is left. The marker may be
        #mid-line (e.g. a trailing // comment); only the rest of the line goes.
        while True:
            start=ln.find('#include')
            if start==-1:
                start=ln.find('//')
            if start==-1:
                start=ln.find('#define')
            if start==-1:
                break
            end=ln.find('\n', start)
            if end==-1: #on the last line without a trailing newline
                ln=ln[0:start]
            else:
                ln=ln[0:start]+ln[end+1:]

        while True: #Remove all C style comments.
            start=ln.find('/*')
            end=ln.find('*/', start+2) #start+2: "/*/" does not close itself
            if start==-1 or end==-1:
                break
            ln=ln[0:start]+ln[end+2:]

        ln=ln.replace('\n',' ')
        ln=ln.replace(';','; ')
        ln=re.sub(r'format,[ ]*\.\.\.','format',ln) #drop variadic args after format
        ln=re.sub(r'[)][ ]*[(]',')(',ln) #replace space between ) and (
        ln=ln.replace('*','* ') #add space after *
        ln=re.sub(r'[ ]*[*]','*', ln) #remove space before *
        ln=re.sub(r'[ ]+',' ', ln) #remove double space
        ln=ln.replace('* *','**') #remove space between **
        lines.append(ln)
    return lines
#parse C structs and its fields, and return variable name and type
def parse_structs(srcdir, files):
    """Parse struct definitions in C headers.

    Parameters:
        srcdir, files: as in parse_file().

    Returns:
        A dict mapping struct name to a dict of {field name: field type}.
        The name is the struct tag ("struct foo_t {...}") or, for an
        anonymous "typedef struct {...} foo_t;", the typedef name. Field
        types have 'struct' and 'const' removed and '*' attached
        (e.g. "dmat*"); only the first word of the type is kept, so multi
        word types such as "unsigned int" are recorded as "unsigned".
        Fields declared as "int a, b;" record only the last name.

        Parsing of a file stops at the first struct with a nested '{'
        (union or inner struct); see the module docstring. A 'struct'
        keyword that is not followed by its own body (forward declaration)
        picks up the next body, giving an entry with a garbage name; these
        entries are harmless because nothing looks them up.

    Raises:
        OSError: if a file cannot be read.
    """
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
                break
            name=ln[struct+1:start]
            end0=ln.find('{', start+1)
            end=ln.find('}', start+1)
            if end==-1: #unterminated body
                break
            if end0!=-1 and end0<end:
                #print(struct, start, end0, end)
                #print (ln[struct:end])
                #print ("embeded {} are not supported")
                break
            fields=ln[start+1:end].split(';')
            if len(name)==0: #typedef struct {...} name;
                end0=ln.find(';', end+1)
                name=ln[end+1:end0]
            name=name.lstrip().rstrip()
            structs[name]=dict()
            for fi in fields:
                tmp=fi.replace('struct','').replace('const','').lstrip().rstrip().split(' ')
                if len(tmp)>=2:
                    structs[name][tmp[-1]]=tmp[0]
    return structs

#parse C functions
def parse_func(srcdir, files):
    """Parse function declarations in C headers.

    Every "name(...)" in the normalized text is examined. The two words
    before '(' are taken as the return type and name, and the text inside
    () as the comma separated arguments. Anything with nested parentheses
    (function pointers, macros with arguments) is skipped, as is anything
    with fewer than two words before '(' or with an argument that is a
    single word other than a number or "..." (this includes "(void)" and
    "()", so functions without arguments are never returned).

    'const' is removed everywhere, including inside identifiers, and
    'const_' type prefixes (const_anyarray) are dropped.

    Parameters:
        srcdir, files: as in parse_file().

    Returns:
        A dict mapping the exported name to [return type, args, C name].
        args is a list of [type, name] pairs; a literal number argument is
        stored as ['number', value]. In mex/aolib.h a declaration
        "type mexname=cname(args);" gives the exported name mexname and the
        C name cname; otherwise both are the same. A later declaration with
        the same exported name replaces an earlier one.

    Raises:
        OSError: if a file cannot be read.
        Exception: if an argument has more than two words after
            normalization (i.e. cannot be split into type and name).
    """
    funcs=dict()
    lines=parse_file(srcdir, files)
    for fn, ln in zip(files, lines):
        ln=ln.replace('const_','').replace('const','')
        ln=ln.replace('* *','**') #remove space between **
        end=0
        while True:
            funname=''
            openp=ln.find('(', end)
            if openp==-1:
                break;
            closep=ln.find(')', openp+1);
            if closep==-1: #unbalanced (. Nothing more to parse.
                break
            openp2=ln.find('(', openp+1);
            if openp2!=-1 and openp2<closep: #nested (). Skip
                closep2=ln.find(')', closep+1)
                if closep2==-1:
                    break
                end=closep2+1
                continue
            elif openp2==closep+1: #()(), e.g. function pointer type. Skip
                closep2=ln.find(')', openp2+2)
                if closep2==-1:
                    break
                end=closep2+1
                continue

            args=ln[openp+1:closep].split(',')
            #Text between the previous match and '(': "... type name"
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
                        except ValueError:
                            msg='ill argument:('+arg+')'
                            funname=''
                elif len(argpair)==2 :
                    funarg.append(argpair)
                else:#failed to parse function
                    msg=fn+': ill argument:('+arg+') in '+ln[end:closep+1].strip()
                    raise(Exception(msg))
            if len(funname)>0:
                funcs[funname]=[funtype, funarg, funname2]
            #else:
                #print("Skipped:", ln[end:closep], msg)
            end=closep+1;
            if end<len(ln) and (ln[end]==';' or ln[end]==' '):
                end+=1
    return funcs
            
