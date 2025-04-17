#Julia wrapper of readbin to read bin files
module Readbin
export readbin, readmany
import SparseArrays
import Glob
import NaturalSort #for natural sorting
import RecursiveArrayTools #convert array of array to higher dimensional array

"C struct for numeric and object arrays"
struct cell
    id::UInt32; id2::UInt32; p::Ptr{Ptr{Cvoid}}; nx::Int; ny::Int; keywords::Cstring;
end
"C struct for sparse arrays"
struct anysp
    id::UInt32; id2::UInt32; px::Ptr{Cdouble}; nx::Int; ny::Int; keywords::Cstring; 
    fp::Ptr{Cvoid}; nzmax::Int; pp::Ptr{Int}; pi::Ptr{Int}
end
"Mapping C ID to Julia type for dense arrays"
numtype=Dict{UInt32, DataType}(
    0x6402=>Float64, 
    0x6403=>Int64, 
    0x6404=>ComplexF64, 
    0x6405=>Int32, 
    0x6408=>Float32,
    0x6409=>ComplexF32,
    0x640A=>Int8,
    0x640B=>Int16,
)
"Mapping C ID to Julia type for sparse arrays"
sparsetype=Dict{UInt32, Tuple{DataType, DataType}}(
    0x6400=>(ComplexF64, Int64),
    0x6401=>(Float64, Int64),
    0x6406=>(ComplexF64, Int32),
    0x6407=>(Float64, Int32),
    0x6432=>(ComplexF32, Int64),
    0x6430=>(Float32, Int64),
    0x6433=>(ComplexF32, Int32),
    0x6431=>(Float32, Int32),
)
"Convert C struct pointer to Julia data. Data is copied to Julia and freed in C side."
function cell2array(x::Ptr{cell})::Any
    if x==C_NULL
        return []
    end
    px=unsafe_load(x,1) #load C struct
    id=px.id&0x64FF
    if haskey(numtype, id) #For numeric array, px.p is pointer to number array
        res=Array{numtype[id]}(undef, (px.nx, px.ny))
        copyto!(res, unsafe_wrap(Array, convert(Ptr{numtype[id]}, px.p), (px.nx, px.ny))) #own is false by default
    elseif haskey(sparsetype, id)
        vtype=sparsetype[id][1] #value type
        itype=sparsetype[id][2] #integer type
        spx=unsafe_load(convert(Ptr{anysp}, x), 1)
        colptr=Vector{itype}(undef, spx.ny+1)
        copyto!(colptr, unsafe_wrap(Array, convert(Ptr{itype}, spx.pp), spx.ny+1))
        colptr.+=1 #convert to 1 based index
        rowval=Vector{itype}(undef, spx.nzmax)
        copyto!(rowval, unsafe_wrap(Array, convert(Ptr{itype}, spx.pi), spx.nzmax))
        rowval.+=1 #convert to 1 based index
        nzval=Vector{vtype}(undef, spx.nzmax)
        copyto!(nzval, unsafe_wrap(Array, convert(Ptr{vtype}, spx.px), spx.nzmax))
        #Create sparse matrix with existing index
        res=SparseArrays.SparseMatrixCSC{vtype,itype}(spx.nx, spx.ny, colptr, rowval, nzval)
    elseif id==0x6421 #x is pointer to a cell, px.p is pointer to pointer array
        res=Array{Any}(undef, (px.nx, px.ny));
        for iy in 1:px.ny
            for ix in 1:px.nx
                ppi=unsafe_load(px.p, ix+(iy-1)*px.nx); #load pointer
                res[ix,iy]=cell2array(convert(Ptr{cell}, ppi));
            end
        end
    else
        println("unknown id=$(px.id)")
        res=[]
    end 
    
    res
end
const aolib=if haskey(ENV, "MAOS_AOLIB")
        ENV["MAOS_AOLIB"]
    else
        "./aolib"
    end
"Read bin files"
function readbin(fn::String)::Any
    x=@ccall aolib.readbin(fn::Cstring)::Ptr{cell};
    res=cell2array(x);
    @ccall aolib.cellfree_do(x::Ptr{cell})::Cvoid;
    res
end
"Read many bin files and return the result in a big array"
function readmany(fns::String)::Any
    res=[]
    fnall=sort(Glob.glob(fns),lt=NaturalSort.natural)
    for fn in fnall
        append!(res, [readbin(fn)])
    end
    try
        res=convert(Array, RecursiveArrayTools.VectorOfArray(res))
    catch e
        println("readmany: unable to convert")
    end
    res,fnall
end
end #end module
