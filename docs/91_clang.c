/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \page page91_clang C Fundementals

   The user is required to have a good understanding of the C fundementals, like
   the pointer, to study and improve the code. Here we list a few to remind the
   user. 

   \section sect-array Arrays
   
   If you declare an array: <code>double a[1024];</code> the symbol actually
   contains the address of a 1024x8 byte memory region that is automatically
   allocated in the stack (malloc, calloc allocations memory in the
   heap). Suppost the address of the memory region is 0xFF00. Assignment
   <code>a[0]=1;</code> puts 1 in memory started at 0xFF00 <code>a[1]=1;</code>
   puts 1 in memory started at 0xFF08.

   Stack is the memory region allocated dynamically when calling a function. All
   variables reside in stack are automatically deleted when the function
   returns. 
   
   Heap is the memory reservoir. The functions malloc, calloc, realloc reserves
   memory in the heap while free removes the reservation. 

   \section sect-pointer Pointers 

   A pointer is a variable that stores the address of a memory region.
   <code>void *p;</code><p> is a basic pointer that just stores the address and can
   not do any thing about the memory it points to. Its value can be assigned to
   any other non const vectr or from any vector (1d or 2d).

   However, a normal pointer like <code>double *p;</code> also contains
   information about the memory it points to. More specifically, the compiler
   know the length of each element so that it can compute the address offset of
   the elements, just like arrays. 

   For example: <code>double *p;</code> if p contains value of 0xFF00.<p>
   <code>p[0]=1;</code>will store a double 1 into memory region at 0xFF00<p>
   <code>p[1]=1;</code>will store a double 1 into memory region at 0xFF08
   because a double requires 8 bytes.<p> <code>p=p+1;</code> p will now contain
   0xFF08, i.e., points to the next element, instead of 0xFF01 <p>
   <code>p[0]=1;</code> is strictly equivalent to <code>*p=1;</code> *p means
   visit the memory pointed by p instead of p itself.  <code>p[1]=1;</code> is
   strictly equivalent to <code>*(p+1)=1;</code><p>

   
   \subsection sect-1d-pointer 1D Pointers

   The simplist usage of pointer is to create and index an 1-D array:
   <code>double *p;</code>.  

   <code>p=mycalloc(1024,double);</code> Allocate a memory block of
   1024x8 byte and store its address in p.<p> <code>p[0]=1;</code> will assign 1 to
   the first number in the array<p>

   This is equivalent to using the arrays.

   <code>double a[1024];</code> Declare an 1-D array.

   <code>a[0]=1;</code > will assign 1 to the first number in the array.  The
   symbol a is like a constant pointer. We can assign it's value to a constant
   pointer: 

   <code>double *const p=a;</code><p> <code>p[0]=1;</code><p>

   \subsection sect-const Const pointers vs. pointers point to const memory

   Beware about a const pointer and pointer that points to constant memory
   regions: 

   <code>const double *p;</code>p is a pointer that points to a memory
   region that can not be changed. Thus <code>p[0]=1;</code> is illegal.

   <code>double * const p;</code>p is a constant pointer that can not have its
   value changed. Thus <code>p=a</code> is illegal, unless you assign at
   declaration: <code>double *const p=a</code>. However, it is legal to modify
   the values in the memory p points to: <code>p[0]=1;</code>.

   \subsection sect-2d-pointer 2D Pointers 

   A 2D pointer is like a 1D pointer. It stores the address of a memory
   region. However, unlike the 1D pointer, the compiler knows some additional
   infos about the 2D pointer. The compiler knows the length of the first
   dimension so that is can compute the memory location given 2D index like
   a[2][3]; <code>double (*p)[2]=calloc(3*2,sizeof(double))</code> can be used
   in the same way as <code>double a[3][2]</code>. The last index changes
   fastest.
   
   The pointer p in our dmat (or cmat) struct points to the allocated memory. It
   can be cast to a two-d vector to facilitize visiting of the elements. The
   cast is conveniently done by PDMAT(a,pa) for dmat a, which simply declares a
   2-d pointer and assign a->p to it: 
   
   <code>double (pa*)[2]=(void*)a->p</code>. a->p is first cast to void *
   because it is a double * and can not be assigned to double (*)[2]
   directly. The array elements (ix,jy) can be accessed by pa[iy][ix]. Notice
   that ix changes fastest.
 */

