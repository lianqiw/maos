\page page91_basics Fundementals

The user is required to have a good understanding of the C fundementals, like
the pointer, to study and improve the code. Here we list a few to remind the
user. 

\section sect-array Memory Access

If you declare an array: <code>double a[1024];</code> the symbol actually
contains the address of a 1024x8 byte memory region that is automatically
allocated in the stack (malloc, calloc allocations memory in the
heap). Suppost the address of the memory region is 0xFF00. Assignment
<code>a[0]=1;</code> puts 1 in memory started at 0xFF00. And <code>a[1]=1;</code>
puts 1 in memory started at 0xFF08.

Stack is the memory region allocated dynamically when calling a function. All
variables reside in stack automatically disappear when the function returns.

Heap is the memory reservoir that persists during the life time of the
executable. The functions \c malloc, \c calloc, \c realloc reserves memory in
the heap while \c free removes the reservation.

\subsection sect-pointer Pointers 

A pointer is a variable that stores the address of a memory region. For
example, <code>void *p;</code> is a basic pointer that just stores the
address and can not do any thing about the memory it points to. Its value can
be assigned to any other non const vector or from any vector (1d or 2d).

However, a normal pointer like <code>double *p;</code> also contains
information about the memory it points to. More specifically, the compiler
know the length of each element so that it can compute the address offset of
the elements, just like arrays. For example: 

    double *p; //assume p=0xFF00.
    p[0]=1;    //stores 1 into memory region at 0xFF00
    p[1]=1;    //stores 1 into memory region at 0xFF08 because a double requires 8 bytes.
    p=p+1;     //now p=0xFF08, i.e., points to the next element, instead of 0xFF01 
    p[i]=1;    //is strictly equivalent to <code>*(p+i)=1;</code> 

\subsection sect-1d-pointer Array Pointers

A pointer can be used to create and index an 1-dimensional (1-d) array in the heap:

    double *p; //decleares a double array pointer.
    p=calloc(1024,sizeof(double)); //Allocate a memory block of 1024x8 byte and store its address in p.
    p[0]=1;    //will assign 1 to the first number in the array.
    
\subsection sect-const Constness

Beware about a const pointer and pointer that points to constant memory
regions:

    const double *p2;  //p2 is a pointer that points to a memory region that can not be changed. 
    p2=p;              //ok
    p2[0]=1;           //illegal.
    double *const p3=p;//p2 is a constant pointer that can not have its value changed. 
    p3=p;              //illegal
    p3[0]=1;           //ok

\subsection sect-2d-pointer Two Dimension Pointers 

A two-dimensional (2-d) pointer is like an array pointer. It stores the
address of a memory region. However, unlike the array pointer, the compiler
knows the length of the first dimension so that is can compute the memory location
given 2-d index like \c a[2][3]; 

    double (*p)[3]=calloc(3*2,sizeof(double)); //can be used in the same way as double a[2][3];

Notice that in 2-d indexing \c a[2][3], the last index changes the fastest.   

\section sect-maos-data MAOS Data Types

\subsection sect-maos-math Fundamental Math Types

MAOS defines and uses a few structs to manage memory for basic math
types. Those types are defined in the \c math folder.

For example, #dmat is used for 1-d or 2-d double array, and #cmat is
used for complex array. There is a predefined macro #P() to help indexing the
arrays:

    dmat *a=dnew(10,4); //allocates memory like double a[4][10]
    P(a,8,3)=1.;        //visits memory like a[3][8]
    P(a,39)=1.;         //visits memory as if it is 1d array.
    dfree(a);           //frees the memory
    cmat *a2=cnew(10,4);
    P(a2,7,3)=1.+1*I;
    cfree(a2);

Sparse arrays are managed using #dsp for double and #csp for complex.

    dsp *a3=dspnew(10,8,50); //10x8 array with maximum of 50 elements.

In addition, #loc_t is used to describe x and y coordinate of uniformly
spaced points. #map_t is derived from #dmat and used to describe 2-d array
of opd or amplitude data with origin and resolution information.

Numerical arrays can be stores in #cell arrays like in matlab. Different
types have their specific cell definitions. For example, #dcell is used to
contain #dmat pointers, and #ccell is used to contain #cmat pointers. All
cell arrays share the same struct layout and can be allocated using
cellnew() and deallocated using cellfree(). However, access elements is
better handled with specific cell types to avoid casting. The #cell array can
contain arbitrary points.

    dcell *c=cellnew(4,4);//allocates a cell array.
    P(c,2,1)=a;           //stores dmat pointer in a cell.
    cellfree(c);

    ccell *c2=cellnew(4,4);
    P(c2,2,3)=a2;
    
    dspcell*c3=cellnew(5,5);
    P(c3,3,2)=a3;
    
    cell *c4=cellnew(5,5);
    P(c4,0,0)=a;
    P(c4,1,0)=a2;
    P(c4,2,0)=a3;


Functions for those fundamental math types are defined using macros that works similarly to the C++ template.

\subsection sect-specific-type Specific Types

Specific calculation related types are often encapsulated in corresponding data types. Those types are defined in the \c lib folder. 

\section sect-bin   Bin file format

The \c bin file format is a simple custom binary format developed specifically to save telemetry data for \c maos. 
It can be converted to \c fits file using the supplied binary \c bin/bin2fits. 
The data in \c bin file are stored with the native endianness to facilitate memory mapping. 

Data in \c bin file are stored as one or multiple consecutive \c blocks :
\verbatim
block1 :self described data
[block2] :(optional) self described data
...
\endverbatim

Each \c block is composed of a \c header and the \c data. 

The \c header contains an optional part followed by a mandatory part. The optional part describes the data:

\verbatim
0x6500: 4 byte uint
length: 8 byte uint equal to the length of the follwoing string
string: variable length, padded to multiple of 8.
length2: 8 byte uint equal to length
0x6500: 4 byte uint
\endverbatim

The mendatory part contains the type and size of \c data
\verbatim
0x6600: 4 byte uint. Dummy, to facilitate memory alignment for mmap.
magic:  4 byte uint. See sys/bin.h for all supported types
nx:     8 byte uint. Inner dimension (fast changing index)
ny:     8 byte uint. Outer dimension (slow changing index)
\endverbatim

The \c data can be either array of fundamental data type or array of other arrays. 
For array of the fundamental data type, it is simply a copy of the data in the native endianness. 
For arrays of other arrays, the \c magic is set to 0x6420 and the data is composed of `nx*ny` \c blocks. 
The \c data part is skipped if either \c nx or \c ny equals to zero.

The \c bin file is readable in Matlab or Python using the supplied \c read or \c readbin routines. In \c maos, it can be either read or mapped to memory. 
For example The atmospheric screens are mapped to memory which enables using extra large atmospheric screens that may not fit in the memory 
and also enables sharing between concurrent maos runs that uses the same seed. 

\section sect-guide Guidelines

When modifying the code, please adhere to the following guidelines.

Memory allocation:

- When declaring a struct, always initialize it to zero by calling calloc or
memset, unless you have a good reason not to.

- Use calloc instead of malloc to initialize the memory unless for large
arrays that you will initialize immediately.

Functions:
    
- The input parameters that is also the output should be grouped in the
beginning, except when there is a main object for the algorithm.

- In any function, there should usually be at most 1 return statement unless
NULL or error code are returned before matural.

- One utility function should be handle only one major mask to maximize
reusability.

- Input arguments to utility functions (in lib folder) should be basic types
to maximize resuability. Input arguments to simulation functions can be
wrapped (like simu) to maximum readability.

Others:

- Avoid non constant static variables/pointers unless you have a good reason
to do so. static variables are bad for multi-threading and hard to free
after usage.

- whenever modify something temporarily for debugging, make a warning. let
it easily identifiable. not make a hidden bug.

- Do not hard code adjustable parameters.

- Always declare unchanged variables to constants. declare pointers as
restrict if possible.

- Do not include system headers in header files unless necessary. include
thoese headers in .c file. 

- Avoid function casting. It will hide data type check and hide bugs.

