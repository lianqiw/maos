function writebin(data,varargin)
%WRITEBIN Write Matlab data to a MAOS .bin file.
%
%   WRITEBIN(DATA, FN) writes DATA to file FN.
%   WRITEBIN(DATA, HEADER, FN) also writes the text HEADER with the data.
%
%   Inputs:
%     DATA    Numeric, logical, or char array (at most 2-d), real or complex,
%             dense or sparse, or a (nested) cell array of these. Real integer
%             arrays are written as integers: uint8, int16, int32, and int64
%             as is (M_UINT8, M_INT16, M_INT32, M_INT64), int8 as int16, uint16
%             and uint32 as the next larger signed type, and uint64 as int64
%             (values must fit). Other dense data is written as double (M_DBL) or
%             double complex (M_CMP). Sparse data is written as
%             M_DSP32/M_CSP32 (32 bit indices) when possible, otherwise as
%             M_DSP64/M_CSP64. An all-zero sparse matrix is written as an
%             empty (0x0) array, as MAOS does. Cells are written as MCC_ANY.
%     HEADER  Text written before the data, normally key=value entries
%             separated by newlines, which MAOS reads as keywords (e.g.
%             sprintf('dx=%g\nox=%g', 0.1, -5); see save_surf.m). For cell
%             DATA, HEADER may be:
%               a char array: header of the whole cell;
%               a cell with numel(DATA) entries: header of each element;
%               a cell with numel(DATA)+1 entries: header of each element,
%               and the last one is the header of the whole cell.
%             Note that READBIN returns either the element headers or, if
%             there are none, the header of the whole cell, not both.
%     FN      Output file name. If it ends with .bin.gz (or .bin.zst) the file
%             is compressed with gzip (or the zstd command). .bin is appended
%             if FN does not end with .bin (a .gz suffix not preceded by .bin is
%             dropped, without compression). The .fits format is not supported.
%
%   Format: see READBIN and sys/bin.c. Each array is written as
%   [M_COMMENT block] uint32 magic, uint64 nx, uint64 ny, data, in native
%   (little endian) byte order.
%
%   Example:
%     writebin(rand(10), 'dx=0.1', 'opd.bin');
%     writebin({1, [2 3]; sparse(eye(3)), 1i}, 'test.bin.gz');
%
%   See also READBIN.
    if nargin==3
        header=varargin{1};
        fn=varargin{2};
    elseif nargin==2
        fn=varargin{1};
        header='';
    else
        error('Usage: writebin(data,header,filename) or writebin(data, filename)');
    end
    if ~ischar(fn) || isempty(fn)
        error('writebin: file name must be a non-empty char array');
    end
    if ~isempty(header) && ~ischar(header) && ~iscell(header)
        error('writebin: header must be a char array or a cell array of char arrays');
    end
    do_gzip=0;%default don't do compress
              %if fn ends with .bin, we don't do gzip
    do_zstd=0;
    if length(fn)>3 &&strcmp(fn(end-3:end),'.bin')
        do_gzip=0;
    end
    %if fn ends with .gz, strip it off.
    if length(fn)>2 && strcmp(fn(end-2:end),'.gz')
        fn=fn(1:end-3);
        do_gzip=1;
    end
    %if fn ends with .zst, strip it off.
    if length(fn)>3 && strcmp(fn(end-3:end),'.zst')
        fn=fn(1:end-4);
        do_zstd=1;
    end
    if length(fn)>4 && strcmp(fn(end-4:end),'.fits')
        warning('writebin:fits', 'writebin: .fits is not supported; writing %s.bin in .bin format', fn);
    end
    %if fn doesn't end with .bin, add .bin
    if ~(length(fn)>3 &&strcmp(fn(end-3:end),'.bin'))
        do_gzip=0;
        do_zstd=0;
        fn=[fn '.bin'];
    end
    fid=fopen(fn,'wb','l');
    if fid<0
        error('writebin: unable to open %s for writing', fn);
    end
    cleanup=onCleanup(@() close_if_open(fid)); %close the file if an error happens
    writebin_do(data,header,fid);
    if fclose(fid) %must close (flush) before compressing
        error('writebin: error closing %s', fn);
    end
    if do_gzip
        try
            gzip(fn); %writes fn.gz
            delete(fn);
        catch
            [status, msg]=system(sprintf('gzip -f "%s"', fn));
            if status
                error('writebin: gzip failed for %s: %s', fn, msg);
            end
        end
    elseif do_zstd
        [status, msg]=system(sprintf('zstd -q -f --rm "%s"', fn)); %writes fn.zst
        if status
            error('writebin: zstd failed for %s: %s', fn, msg);
        end
    end
function close_if_open(fid)
    if any(fopen('all')==fid)
        fclose(fid);
    end
function writeheader_do(header, fid)
    M_COMMENT=25856;
    if ~ischar(header)
        error('writebin: header must be char type (a cell header is only allowed for cell data)');
    end
    fwrite(fid,uint32(M_COMMENT),'uint32');
    fwrite(fid,uint64(length(header)+1),'uint64');
    fwrite(fid,header,'char');
    fwrite(fid,0,'char');
    fwrite(fid,uint64(length(header)+1),'uint64');
    fwrite(fid,uint32(M_COMMENT),'uint32');
    
function writebin_do(data,header,fid)
    M_CSP64=25600;
    M_SP64=25601;
    M_CSP32=25606;
    M_SP32=25607;
    M_DBL=25602;
    M_INT64=25603;
    M_CMP=25604;
    M_INT32=25605;
    M_UINT8=25610;
    M_INT16=25611;
    MCC_ANY=25633;
    MAT_SP=65281;
    MAT_CSP=65282;
    M_COMMENT=25856;
    if length(size(data))>2
        error('writebin: We do not handle more than 2 dimensions');
    end
    if ~iscell(data) && ~isnumeric(data) && ~islogical(data) && ~ischar(data)
        error('writebin: unsupported data type %s', class(data));
    end
    if iscell(data)
        if iscell(header) && ~isempty(header) && numel(header)~=numel(data) && numel(header)~=numel(data)+1
            error('writebin: a cell header must have numel(data) or numel(data)+1 entries; got %d for %d', numel(header), numel(data));
        end
        if ~isempty(header) 
            if ~iscell(header) %header is not cell.
                writeheader_do(header, fid);
            else
                if numel(header)==numel(data)+1
                    writeheader_do(header{numel(data)+1},fid);
                end
            end
        end
        fwrite(fid,MCC_ANY,'uint32');
        fwrite(fid,size(data,1),'uint64');
        fwrite(fid,size(data,2),'uint64');
        if ~isempty(header) && iscell(header)
            for ii=1:numel(data)
                writebin_do(data{ii},header{ii},fid);
            end
        else
            for ii=1:numel(data)
                writebin_do(data{ii},[],fid);
            end
        end
    else
        %write the header.
        if ~isempty(header)
            writeheader_do(header, fid);
        end
        if issparse(data)
            [Ir,Ic,P]=find(data);
            nz=length(P); %was nzmax(data), which may exceed the number of entries written
            P=double(P); %logical sparse
            if isreal(P)
                complex=0;
            else
                complex=1;
                P2=zeros(2,length(P));
                P2(1,:)=real(P);
                P2(2,:)=imag(P);
                P=P2(:);
            end
            if nz<intmax('int32') && size(data,1)<intmax('int32') %use 32bit (MAOS reads as signed int) to save
                use_32=1;
                if complex
                    magic=M_CSP32;
                else
                    magic=M_SP32;
                end
            else
                use_32=0;
                if complex
                    magic=M_CSP64;
                else
                    magic=M_SP64;
                end
            end
            nx=size(data,1);
            ny=size(data,2);
            if nz==0 %as spwritedata in math/spbin.c: an empty sparse matrix is written as 0x0
                nx=0;
                ny=0;
            end
            fwrite(fid,magic,'uint32');
            fwrite(fid,nx,'uint64');
            fwrite(fid,ny,'uint64');
            if nz>0
                fwrite(fid,nz,'uint64');
                Jc=zeros(ny+1,1);
                pos=1;
                for icol=1:ny
                    Jc(icol)=pos;
                    while pos<=nz && Ic(pos)==icol
                        pos=pos+1;
                    end
                end
                Jc(ny+1)=pos;
                %first convert to 0 index
                Ir=Ir-1;
                Jc=Jc-1;
                if use_32
                    Ir=uint32(Ir);
                    Jc=uint32(Jc);
                    fwrite(fid,Jc,'uint32');
                    fwrite(fid,Ir,'uint32');
                else
                    Ir=uint64(Ir);
                    Jc=uint64(Jc);
                    fwrite(fid,Jc,'uint64');
                    fwrite(fid,Ir,'uint64');
                end
                fwrite(fid,P,'double');
            end
        elseif isinteger(data) && isreal(data) %integers of the same or larger size
            switch class(data)
             case 'uint8'
              magic=M_UINT8; prec='uint8';
             case {'int16', 'int8'}
              magic=M_INT16; prec='int16';
             case {'int32', 'uint16'}
              magic=M_INT32; prec='int32';
             otherwise %int64, uint32, uint64
              if isa(data, 'uint64') && any(data(:)>uint64(intmax('int64')))
                  error('writebin: uint64 value does not fit in int64');
              end
              magic=M_INT64; prec='int64';
            end
            fwrite(fid,magic,'uint32');
            fwrite(fid,size(data,1),'uint64');
            fwrite(fid,size(data,2),'uint64');
            fwrite(fid,data,prec);
        elseif isreal(data) %other numeric, logical, and char are converted to double
            fwrite(fid,M_DBL,'uint32');
            fwrite(fid,size(data,1),'uint64');
            fwrite(fid,size(data,2),'uint64');
            fwrite(fid,data,'double');
        else %complex
            fwrite(fid,M_CMP,'uint32');
            fwrite(fid,size(data,1),'uint64');
            fwrite(fid,size(data,2),'uint64');
            tmp=zeros(2,numel(data));
            tmp(1,:)=real(data(:));
            tmp(2,:)=imag(data(:));
            fwrite(fid,tmp,'double');
            clear tmp;
        end 
    end