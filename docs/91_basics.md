\page page91_basics Fundamentals

To study and improve the code you need a good understanding of C
fundamentals, such as pointers. This page reviews a few of them, then
describes the data types, the \c bin file format, and the coding guidelines
used in MAOS.

\section sect-array Memory Access

If you declare an array <code>double a[1024];</code>, the symbol \c a refers to
a 1024x8 byte memory region that is automatically allocated on the stack
(\c malloc and \c calloc allocate memory on the heap). Suppose the address of
the memory region is 0xFF00. The assignment <code>a[0]=1;</code> puts 1 in the
memory starting at 0xFF00, and <code>a[1]=1;</code> puts 1 in the memory
starting at 0xFF08.

The stack is the memory region allocated automatically when a function is
called. All variables on the stack disappear when the function returns.

The heap is the memory reservoir that persists during the lifetime of the
program. The functions \c malloc, \c calloc, and \c realloc reserve memory on
the heap, while \c free releases the reservation.

\subsection sect-pointer Pointers

A pointer is a variable that stores the address of a memory region. For
example, <code>void *p;</code> is a basic pointer that only stores the
address; it cannot be used to access the memory it points to. A <code>void
*</code> pointer can be assigned to or from any other (non-const) data pointer
without a cast in C.

A typed pointer like <code>double *p;</code> also carries information about
the memory it points to. The compiler knows the size of each element, so it can
compute the address offset of the elements, just like arrays. For example:
```
    double *p; //assume p=0xFF00.
    p[0]=1;    //stores 1 into memory at 0xFF00
    p[1]=1;    //stores 1 into memory at 0xFF08 because a double requires 8 bytes.
    p=p+1;     //now p=0xFF08, i.e., points to the next element, instead of 0xFF01
    p[i]=1;    //is strictly equivalent to *(p+i)=1;
```
\subsection sect-1d-pointer Array Pointers

A pointer can be used to create and index a 1-dimensional (1-d) array on the heap:
```
    double *p; //declares a double array pointer.
    p=calloc(1024,sizeof(double)); //Allocate a memory block of 1024x8 bytes and store its address in p.
    p[0]=1;    //assigns 1 to the first number in the array.
    free(p);   //release the memory.
```
\subsection sect-const Constness

Distinguish a pointer to constant memory from a constant pointer:
```
    const double *p2;  //p2 points to memory that cannot be changed through p2.
    p2=p;              //ok
    p2[0]=1;           //illegal.
    double *const p3=p;//p3 is a constant pointer: its value cannot be changed.
    p3=p;              //illegal
    p3[0]=1;           //ok
```
In MAOS, input arguments that are not modified are declared \c const, for
example <code>void dcp(dmat **out, const dmat *in);</code>.

\subsection sect-2d-pointer Two Dimension Pointers

A pointer to an array (a 2-d pointer) also stores the address of a memory
region. Unlike the plain array pointer, the compiler knows the length of the
last dimension, so it can compute the memory location of a 2-d index like
\c a[2][3]:
```
    double (*p)[3]=calloc(3*2,sizeof(double)); //can be used in the same way as double a[2][3];
```
Notice that in C 2-d indexing \c a[2][3], the last index changes the fastest.
The MAOS array types below use the opposite (column-major, Fortran or Matlab)
order: the first index changes the fastest.

\section sect-maos-data MAOS Data Types

\subsection sect-maos-math Fundamental Math Types

MAOS uses a few structs to manage memory for basic math types. They are
defined in \c math/type.h and the functions are declared in the \c math folder,
except #loc_t, #pts_t, and #map_t, which are defined in \c lib/loc.h and
\c lib/map.h.
Each array struct records its type (a magic number, see \ref sect-bin), the
data pointer \c p, the dimensions \c nx (rows, the fast changing index) and
\c ny (columns), optional \c keywords (a string saved with the data), and
reference counting information so that several arrays can share the same
memory.

| Type | Element | Description |
|------|---------|-------------|
| #dmat | \c real | real array (double, or float with `--disable-double`) |
| #cmat | \c comp | complex array (double complex, or float complex with `--disable-double`) |
| #smat | \c float | single precision real array |
| #zmat | \c fcomplex | single precision complex array |
| #lmat | \c long | integer array |
| #dsp, #csp | \c real, \c comp | sparse matrix in compressed sparse column format |
| #ssp, #zsp | \c float, \c fcomplex | single precision sparse matrix |
| #loc_t | \c real | x and y coordinates of a list of points (\c locx, \c locy, \c nloc) |
| #map_t | \c real | 2-d OPD or amplitude map with origin (\c ox, \c oy) and sampling (\c dx, \c dy) |
| #cell | any | cell array of pointers to any of the above |

\c real and \c comp are defined in \c math/numtype.h. They are \c double and
\c dcomplex by default. When configured with `--disable-double`, the macro
\c CPU_SINGLE is set to 1 and \c real and \c comp become \c float and
\c fcomplex. The \c d and \c c prefixed functions then operate on single
precision data, so code that uses \c real works in both modes.

Use the macros defined in \c math/type.h to access the arrays:

| Macro | Meaning |
|-------|---------|
| \c P(A) | the data pointer \c A->p |
| \c P(A,i) | element \c i, treating the array as 1-d |
| \c P(A,ix,iy) | element at row \c ix and column \c iy, i.e., \c A->p[ix+iy*A->nx] |
| \c PR(A,ix,iy) | like \c P(), but an index is wrapped to 0 when that dimension has size 1 |
| \c PCOL(A,iy) | pointer to the first element of column \c iy |
| \c PN(A), \c NX(A), \c NY(A) | number of elements, rows, and columns; 0 if \c A is \c NULL |

When the code is compiled with `--enable-debug`, \c P() checks that the index
is within the array bounds and calls error() otherwise.

```
    dmat *a=dnew(10,4); //allocates memory like double a[4][10]
    P(a,8,3)=1.;        //visits memory like a[3][8]
    P(a,39)=1.;         //visits memory as if it is a 1d array.
    dfree(a);           //frees the memory and sets a to NULL.
    cmat *a2=cnew(10,4);
    P(a2,7,3)=COMPLEX(1,1);
    cfree(a2);
```
Sparse arrays are managed using #dsp for real and #csp for complex values.
```
    dsp *a3=dspnew(10,8,50); //10x8 sparse matrix with room for 50 non-zero elements.
    dspfree(a3);
```
In addition, #loc_t describes the x and y coordinates of a set of points,
usually on a regular grid with spacing \c dx and \c dy. #map_t describes a 2-d
array of OPD or amplitude data together with its origin, sampling, and
conjugation height. Its first part has the same layout as #dmat, so a map can
be passed to #dmat functions via \c map->dmat; likewise \c loc->dmat views the
coordinates as an \c nloc x 2 #dmat. \c math/type.h only declares them, so
that the generic routines such as \c writebin() and \c cellfree() accept them,
and checks at compile time that their layout agrees with #dmat.

Numerical arrays can be stored in #cell arrays like in Matlab. Each type has
its own cell type: for example, #dcell contains #dmat pointers, #ccell
contains #cmat pointers, and #dspcell contains #dsp pointers. Cells of cells
(#dccell, ...) and cells of maps and locs (#mapcell, #loccell) also exist. All
cell arrays share the same struct layout. Allocate them with the typed
constructors, such as \c dcellnew() or \c ccellnew() (macros that call
cellnew() and cast the result), and free them with cellfree() or the typed
versions such as dcellfree(). Freeing a cell also frees its elements. Accessing
elements is best done through the specific cell types to avoid casting. The
generic #cell type can hold arrays of different types.
```
    dcell *c=dcellnew(4,4);//allocates a cell array of dmat.
    P(c,2,1)=dnew(10,4);   //stores a dmat pointer in the cell.
    dcellfree(c);          //frees the cell and the dmat it contains.

    dspcell *c3=dspcellnew(5,5);
    P(c3,3,2)=dspnew(10,8,50);
    dspcellfree(c3);

    cell *c4=cellnew(5,5); //generic cell. Elements must be cast to cell*.
    P(c4,0,0)=(cell*)dnew(10,4);
    P(c4,1,0)=(cell*)cnew(10,4);
    cellfree(c4);
```
Many functions accept any array type through the transparent unions
\c anyarray, \c anycell, and their \c const variants defined in
\c math/type.h, so no cast is needed when calling them, for example
writebin() and cellfree().

\subsection sect-type-generic Type Generic Implementation

Functions for the fundamental math types are written once and compiled for
each type, similar to C++ templates. \c math/defs.h defines the macro \c X(name)
that expands to \c dname, \c sname, \c cname, \c zname, or \c lname, and the
element type \c T and its real counterpart \c R, depending on the flags
\c COMP_COMPLEX, \c COMP_SINGLE, and \c COMP_LONG. The files \c math/dmath.c,
\c smath.c, \c cmath.c, \c zmath.c, and \c lmath.c each set these flags and
include the same template sources (\c mat.c, \c matmath.c, \c blas.c,
\c sp.c, ...). To add a function, write it once in the template file as
\c X(foo) using \c T and \c R, and declare it in the corresponding header
macro (for example \c math/mat.h). Do not edit the per-type wrapper files.

\subsection sect-specific-type Specific Types

Types used by specific calculations are defined together with the code that
uses them in the \c lib folder, for example \c propdata_t for ray tracing
(\c lib/accphi.h), \c locfft_t for PSF computation (\c lib/locfft.h), and
\c hyst_t for DM hysteresis (\c lib/hyst.h). The simulation specific structs
are in \c maos/parms.h and \c maos/types.h (see \ref page93_maos).

\section sect-bin Bin file format

The \c bin file format is a simple binary format developed to save telemetry
data for \c maos. The data are stored with the native endianness to allow
memory mapping. MAOS writes and reads data with writebin() and readbin() (or
typed versions such as dread() and dcellread()). The format is selected by the
file name suffix (see \c zfopen() in \c sys/bin.c):

- \c .bin: the \c bin format described below. If a file name has no recognized
  suffix, \c .bin is appended.
- \c .fits: the FITS format. Arrays are saved as FITS images and the keywords
  as FITS header cards. The elements of a cell array are saved as consecutive
  HDUs (the primary image followed by extensions).
- \c .bin.gz, \c .fits.gz: gzip compressed. \c .bin.zst, \c .fits.zst: zstd
  compressed (when compiled with zstd support). When reading, compression is
  detected from the file content (magic bytes), not the name.

The supplied program \c bin/bin2fits converts \c bin files to \c fits files.

Data in a \c bin file are stored as one or multiple consecutive \c blocks:
```
block1 :self described data
[block2] :(optional) self described data
...
```

Each \c block is composed of a \c header and the \c data.

The \c header contains an optional part followed by a mandatory part. The
optional part contains the keywords string that describes the data:

```
0x6500: 4 byte uint (M_COMMENT)
length: 8 byte uint equal to the length of the following string
string: variable length, null terminated and padded with zeros to a multiple of 8 bytes.
length: 8 byte uint, same as above
0x6500: 4 byte uint (M_COMMENT)
```

The string consists of \c key=value entries separated by new lines, just like
the configuration files. If several such parts are present, the strings are
concatenated when read.

The mandatory part contains the type and size of the \c data:
```
0x6600: 4 byte uint (M_SKIP). Padding to facilitate memory alignment for mmap.
magic:  4 byte uint. See sys/bin.h for all supported types
nx:     8 byte uint. Inner dimension (fast changing index)
ny:     8 byte uint. Outer dimension (slow changing index)
```

Commonly used \c magic numbers (enum in \c sys/bin.h):

| %Magic | Name | Data |
|-------|------|------|
| 0x6402 | \c M_DBL | double |
| 0x6408 | \c M_FLT | float |
| 0x6404 | \c M_CMP | double complex |
| 0x6409 | \c M_ZMP | float complex |
| 0x6403, 0x6405 | \c M_INT64, \c M_INT32 | 64 or 32 bit integers |
| 0x640A, 0x640B | \c M_INT8, \c M_INT16 | 8 or 16 bit integers |
| 0x6401, 0x6400 | \c M_DSP64, \c M_CSP64 | sparse double/complex with 64 bit index |
| 0x6407, 0x6406 | \c M_DSP32, \c M_CSP32 | sparse double/complex with 32 bit index |
| 0x6430--0x6433 | \c M_SSP64, \c M_SSP32, \c M_ZSP64, \c M_ZSP32 | single precision sparse |
| 0x6421 | \c MCC_ANY | cell array |
| 0x64FF | \c M_EOD | end of data marker, used in socket transfers |

The \c data can be either an array of a fundamental data type or an array of
other arrays:

- For an array of a fundamental data type, it is simply a copy of the
  \c nx*ny elements in the native endianness.
- For a sparse matrix, it contains the number of non-zero elements \c nzmax
  (8 bytes), the column pointers (\c ny+1 integers), the row indices (\c nzmax
  integers), and the values (\c nzmax elements).
- For an array of other arrays (a cell), the \c magic is set to 0x6421
  (\c MCC_ANY) and the data is composed of \c nx*ny \c blocks. An empty
  element is saved as a block with zero dimensions.

The \c data part is empty if either \c nx or \c ny equals zero.

The \c bin file is readable in Python with \c scripts/readbin.py, in Matlab
with \c scripts/readbin.m or the \c read mex routine built from \c mex/read.c,
in Julia with \c scripts/readbin.jl, and in IDL with \c scripts/readbin.pro. In
\c maos, a file can be either read or mapped into memory. For example, the
atmospheric screens can be mapped into memory, which enables using very large
screens that may not fit in memory and sharing between concurrent maos runs
that use the same seed.

\section sect-guide Guidelines

When modifying the code, please adhere to the following guidelines.

Memory allocation:

- When declaring a struct, always initialize it to zero by calling \c calloc
  (or \c mycalloc) or \c memset, unless you have a good reason not to.

- Use \c calloc instead of \c malloc to initialize the memory, except for large
  arrays that you will initialize immediately.

- Free everything you allocate. Arrays are freed with the matching function
  (\c dfree, \c cfree, \c cellfree, ...), which also sets the pointer to
  \c NULL. Memory accounting (see \ref page95_profiling) reports memory that
  is not freed.

Functions:

- Input parameters that are also output should be grouped at the beginning,
  except when there is a main object for the algorithm.

- A function should usually have at most one return statement, except when
  returning \c NULL or an error code early.

- A utility function should handle only one major task to maximize
  reusability.

- Input arguments to utility functions (in the \c lib folder) should be basic
  types to maximize reusability. Input arguments to simulation functions (in
  the \c maos folder) can be wrapped (like \c parms_t or \c sim_t) to improve
  readability.

Others:

- Avoid non-constant static variables/pointers unless you have a good reason
  to do so. Static variables are bad for multi-threading and hard to free
  after use.

- Whenever you modify something temporarily for debugging, emit a
  \c warning() so that it is easily identifiable and does not become a hidden
  bug.

- Do not hard code adjustable parameters. Make them configuration parameters
  (see \ref page93_maos).

- Always declare unchanged variables as \c const. Declare pointers as
  \c restrict if possible.

- Do not include system headers in header files unless necessary. Include
  those headers in the \c .c file.

- Avoid function pointer casting. It hides data type checks and hides bugs.

- Use the logging macros (\c info, \c warning, \c error, \c dbg, ...) instead
  of \c printf so that the output respects \c MAOS_LOG_LEVEL and goes to the
  log file.

- The code must also compile as C++ (CI builds it with \c CC=clang++), so do
  not use C++ keywords such as \c new or \c class as identifiers and cast
  explicitly where C++ requires it.

- Each source file carries the GPL copyright header (see \c scripts/addcopy.sh).
