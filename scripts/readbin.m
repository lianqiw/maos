function [res header]=readbin(fn)
%READBIN Read a MAOS .bin or .fits file into Matlab.
%
%   RES=READBIN(FN) reads the file FN written by MAOS (sys/bin.c) or by
%   WRITEBIN and returns the data.
%   [RES, HEADER]=READBIN(FN) also returns the text header (keywords) stored
%   with the data.
%
%   Input:
%     FN  File name. It may end with .bin, .fits, .bin.gz, .fits.gz, .bin.zst
%         or .fits.zst. If FN does not exist, the suffixes are appended and
%         tried in that order (like MAOS). The format (.bin or .fits) is
%         determined from the name; gzip or zstd compression is detected from
%         the file content and the file is decompressed to a temporary file
%         (.zst requires the zstd command). The original file is not modified.
%
%   Outputs:
%     RES     Dense arrays are returned as double, single, int64, int32,
%             int16, or uint8 according to the stored type (real or complex).
%             FITS BSCALE/BZERO are not applied.
%             Sparse arrays are returned as Matlab (double) sparse matrices.
%             Cell arrays are returned as nested Matlab cell arrays. map_t and
%             loc_t arrays are returned as plain matrices; their parameters
%             are in HEADER. Empty arrays are returned as [].
%             A .fits file with more than one HDU is returned as an Nx1 cell
%             with one entry per HDU.
%     HEADER  The header string (key=value entries separated by newlines).
%             For a cell array, if any element has a header, HEADER is a cell
%             of the same size with the header of each element (the header of
%             the cell itself is then dropped); otherwise it is the header of
%             the cell. For a multi-HDU .fits file it is a cell with one entry
%             per HDU. For .fits files, all cards other than the structural
%             ones are returned, one per line: COMMENT cards as their text,
%             CONTINUE cards appended to the previous line, and other cards
%             as written (e.g. 'DX      =   0.0625').
%
%   File format (.bin, little endian, see sys/bin.c and sys/bin.h): each
%   array is optionally preceded by headers
%     uint32 M_COMMENT(0x6500), uint64 n, char[n] text, uint64 n, uint32 M_COMMENT
%   then an optional uint32 M_SKIP(0x6600) padding, then
%     uint32 magic, uint64 nx, uint64 ny, data.
%   Dense data is nx*ny elements in column major order (complex as
%   interleaved real/imag). Sparse data (if nx*ny>0) is uint64 nz, column
%   pointers Jc(ny+1), row indices Ir(nz) (0 based, 32 or 64 bit), and values.
%   A cell (MCC_ANY 0x6421) is followed by its nx*ny elements.
%
%   Example:
%     [psf, hdr]=readbin('evlpsfcl_1_x0_y0.fits');
%     res=readbin('Res_1');   %tries Res_1.bin, Res_1.fits, ...
%
%   See also WRITEBIN.

    if ~ischar(fn) || isempty(fn)
        error('readbin: file name must be a non-empty char array');
    end
    if strncmp(fn, '~/', 2) %expand home directory as MAOS does
        fn=[getenv('HOME') fn(2:end)];
    end
    if exist(fn, 'file')~=2 %not found or a directory
        %same order as zfopen in sys/bin.c
        suf={'.bin', '.fits', '.bin.gz', '.fits.gz', '.bin.zst', '.fits.zst'};
        found=0;
        for is=1:length(suf)
            if exist([fn, suf{is}], 'file')==2
                fn=[fn, suf{is}];
                found=1;
                break;
            end
        end
        if ~found
            error('readbin: %s not found (also tried .bin, .fits, .bin.gz, .fits.gz, .bin.zst, .fits.zst)', fn);
        end
    end
    isfits=istype(fn, '.fits') || istype(fn, '.fits.gz') || istype(fn, '.fits.zst');

    %MATLAB does not support transparent decompression with fopen. 
    %So decompress the file to a temporary file first. Compression is
    %detected from the content, like zfopen in sys/bin.c.
    fn0=fn;
    [fn, istmp]=decompress(fn);
    if isfits
        fid=fopen(fn,'rb','b');%handle big endian
    else
        fid=fopen(fn,'rb','l');%MAOS writes native (little endian) byte order
    end
    if fid<0
        if istmp
            delete_tmp(fn);
        end
        error('readbin: unable to open %s for reading', fn0);
    end
    %close the file (and delete the temporary file) on return or error
    cleanup=onCleanup(@() close_file(fid, fn, istmp));
    if isfits %we assembly multiple arrays to a cell
        count=0;
        err=0;
        while ~feof(fid) && ~err
            count=count+1;
            [res{count,1} header{count,1} err]=readbin_do(fid, isfits);
        end
        %last cell is fake when err happens
        if err
            res=res(1:end-1);
            header=header(1:end-1);
        end
        if length(res)==1
            res=res{1};
            header=header{1};
        end
    else
        [res header err]=readbin_do(fid, isfits);
        if err
            warning('readbin:read', 'readbin: %s is empty or truncated', fn0);
        end
    end
    %fid and the temporary file are closed/deleted by onCleanup.
function [fn2, istmp]=decompress(fn)
%Decompress fn to a temporary file if it is gzip or zstd compressed.
%Returns fn unchanged (istmp=0) otherwise.
    fn2=fn;
    istmp=0;
    fid=fopen(fn, 'rb');
    if fid<0
        error('readbin: unable to open %s for reading', fn);
    end
    magic=fread(fid, 4, 'uint8=>uint8')';
    fclose(fid);
    if length(magic)>=2 && magic(1)==31 && magic(2)==139 %gzip: 1f 8b
        tmpdir=tempname;
        mkdir(tmpdir);
        try
            out=gunzip(fn, tmpdir);
        catch err
            rmdir(tmpdir, 's');
            error('readbin: unable to gunzip %s: %s', fn, err.message);
        end
        fn2=out{1};
        istmp=1;
    elseif length(magic)>=4 && all(magic==uint8([40 181 47 253])) %zstd: 28 b5 2f fd
        tmpdir=tempname;
        mkdir(tmpdir);
        fn2=fullfile(tmpdir, 'readbin.tmp');
        [status, msg]=system(sprintf('zstd -q -d -f -o "%s" "%s"', fn2, fn));
        if status
            rmdir(tmpdir, 's');
            error('readbin: unable to decompress %s with zstd: %s', fn, msg);
        end
        istmp=1;
    end
function close_file(fid, fn, istmp)
    fclose(fid);
    if istmp
        delete_tmp(fn);
    end
function delete_tmp(fn2)
%Delete the temporary decompressed file and its directory.
    tmpdir=fileparts(fn2);
    if exist(fn2, 'file')
        delete(fn2);
    end
    if exist(tmpdir, 'dir')
        rmdir(tmpdir, 's');
    end
function res=istype(fn, suffix)
    if length(fn)>length(suffix) && strcmp(fn(end-length(suffix)+1:end), suffix)
        res=1;
    else
        res=0;
    end
function [magic, nx, ny, header]=readfits_header(fid)
%Read the header of one FITS HDU (as read_fits_header in sys/bin.c).
%Every header line (card) is 80 bytes. Total 36 lines or 2880 bytes per page.
%Returns magic=0 at end of file. Complex data is stored by MAOS with NAXIS=3
%and NAXIS1=2 (interleaved real/imag).
    END=0;
    page=0;
    bitpix=0;
    nc=1; %2 for complex
    nx=0;
    ny=0;
    header='';
    while ~END
        start=1;
        if page==0
            res=fread(fid, 80, 'char*1=>char*1')';
            if length(res)<80
                END=1;
                break;
            end
            if ~strncmp(res, 'SIMPLE', 6) && ~strncmp(res, 'XTENSION= ''IMAGE', 16)
                warning('readbin:fits', 'readbin: garbage in fits file; ignore the rest');
                END=1;
                break;
            end
            bitpix=read_card_num(fid);
            naxis=read_card_num(fid);
            if naxis>3
                error('readbin: fits files with NAXIS=%d are not supported', naxis);
            end
            if naxis>2
                nc=read_card_num(fid);
                if nc~=1 && nc~=2
                    error('readbin: fits files with NAXIS=3 and NAXIS1=%d are not supported', nc);
                end
            end
            if naxis>0
                nx=read_card_num(fid);
                ny=1; %1-d data
            end
            if naxis>1
                ny=read_card_num(fid);
            end
            start=start+3+naxis;
        end
        for i=start:36
            res=fread(fid, 80, 'char*1=>char*1')';
            if length(res)<80
                error('readbin: fits file is truncated');
            end
            if strncmp(res, 'END', 3)
                END=1;
            elseif ~END && ~strncmp(res, 'EXTEND', 6) && ~strncmp(res, 'PCOUNT', 6) ...
                    && ~strncmp(res, 'GCOUNT', 6)
                %As in sys/bin.c: keep COMMENT text, append CONTINUE text to
                %the previous line, and keep other cards as is.
                iscont=strncmp(res, 'CONTINUE', 8);
                if strncmp(res, 'COMMENT', 7) || iscont
                    res=res(11:end);
                end
                res=deblank(res);
                if ~isempty(res)
                    if ~isempty(header) && ~iscont
                        header=[header char(10)];
                    end
                    header=[header res];
                end
            end
        end
        page=page+1;
    end
    switch (bitpix)
     case -32
      if nc==2
          magic=25609; %M_ZMP
      else
          magic=25608; %M_FLT
      end
     case -64
      if nc==2
          magic=25604; %M_CMP
      else
          magic=25602; %M_DBL
      end
     case 64
      magic=25603; %M_INT64
     case 32
      magic=25605; %M_INT32
     case 16
      magic=25611; %M_INT16
     case 8
      magic=25610; %M_UINT8
     case 0
      magic=0;
     otherwise
      error('readbin: unsupported fits BITPIX=%d', bitpix);
    end
function val=read_card_num(fid)
%Read one 80 byte FITS card and return its numerical value (columns 11-30).
    res=fread(fid, 80, 'char*1=>char*1')';
    if length(res)<80
        error('readbin: fits file is truncated');
    end
    val=str2double(res(11:min(30,end)));
    if isnan(val)
        error('readbin: unable to parse fits card: %s', res);
    end

function [magic, nx, ny, header]=readbin_header(fid)
%Read the header of one .bin array (as read_bin_header in sys/bin.c): skip
%M_SKIP padding, concatenate M_COMMENT blocks, then read magic, nx, ny.
%Returns magic=0 at end of file or at the end of data marker (M_EOD).
    M_SKIP=26112;
    M_COMMENT=25856;
    M_EOD=25855;

    magic=read_magic(fid, M_SKIP);
    header='';
    while magic==M_COMMENT
        nlen=fread(fid, 1, 'uint64=>double');
        if isempty(nlen)
            error('readbin: header verification failed (file is truncated)');
        end
        str=fread(fid, nlen, 'char*1=>char*1')';
        nlen2=fread(fid, 1, 'uint64=>double');
        magic2=fread(fid, 1, 'uint32=>uint32');
        if length(str)~=nlen || isempty(nlen2) || nlen~=nlen2 || isempty(magic2) || magic2~=M_COMMENT
            error('readbin: header verification failed');
        end
        %MAOS pads the text with nulls to a multiple of 8 bytes. The C reader
        %stops at the first null, so do the same.
        inull=find(str==0, 1);
        if ~isempty(inull)
            str=str(1:inull-1);
        end
        header=[header char(str)];
        magic=read_magic(fid, M_SKIP);
    end
    if isempty(magic) || magic==M_EOD
        magic=0;
        nx=0;
        ny=0;
        return;
    end
    nx=fread(fid,1,'uint64=>double');
    ny=fread(fid,1,'uint64=>double');
    if isempty(nx) || isempty(ny)
        error('readbin: file is truncated in the array header');
    end
    if magic==0 && (nx==0 || ny==0)
        %MAOS writes magic 0 for the empty elements of a cell whose elements
        %are all empty (writedata in math/cell.c). Treat it as empty double.
        magic=uint32(25602);
    end
function magic=read_magic(fid, M_SKIP)
%Read a magic number, skipping any M_SKIP padding. Empty at end of file.
    magic=fread(fid,1,'uint32=>uint32');
    while ~isempty(magic) && bitand(magic, M_SKIP)==M_SKIP
        magic=fread(fid,1,'uint32=>uint32');
    end
function [res, count]=fread_fits(varargin)
%fread for fits: also read and discard the padding to a multiple of 2880 bytes.
    [res, count]=fread(varargin{:});
    fid=varargin{1};
    tmp=whos('res');
    nbyte=tmp.bytes;
    nbyte2=ceil(nbyte/2880)*2880;
    if nbyte<nbyte2 %read in and discard the padding
        tmp=fread(fid, nbyte2-nbyte, 'char*1=>char*1');
    end
function [res header err]=readbin_do(fid, isfits)
%Read one array (recursively for cells). err=1 at end of file/data.
    M_CSP64=25600;
    M_DSP64=25601;
    M_DBL=25602;
    M_INT64=25603;
    M_CMP=25604;
    M_INT32=25605;
    M_CSP32=25606;
    M_DSP32=25607;
    M_FLT=25608;
    M_ZMP=25609;
    M_UINT8=25610;
    M_INT16=25611;
    
    M_SSP64=25648;
    M_SSP32=25649;
    M_ZSP64=25650;
    M_ZSP32=25651;
    
    MC_CSP=25616;
    MC_SP=25617;
    MC_DBL=25618;
    MC_INT64=25619;
    MC_CMP=25620;
    MC_INT32=25621;
    
    MCC_ANY=25633;
    MCC_DBL=25634;
    MCC_CMP=25636;


    MAT_SP=65281;
    MAT_CSP=65282;
    if isfits
        [magic, nx, ny, header]=readfits_header(fid);
        zfread=@fread_fits;
    else
        [magic, nx, ny, header]=readbin_header(fid);
        zfread=@fread;
    end
    if isempty(magic) || magic<=0
        err=1;
        res=[];
        return;
    else
        magic=bitand(magic, 65535);
        err=0;
    end
    if nx==0 || ny==0
        res=[];
        return;
    end
    convert_sparse=0;
    cnt=-1;      %number of elements read by zfread, to detect truncated files
    nread=nx*ny; %number of elements expected
    switch magic
     case {MCC_ANY, MCC_DBL, MCC_CMP, MC_CSP, MC_SP, MC_DBL, MC_CMP, MC_INT32, MC_INT64}
      res=cell(nx,ny);
      header2=cell(nx,ny);
      nheaderstr=0;
      for ii=1:nx*ny
          [res{ii} header2{ii} err]=readbin_do(fid, isfits);
          if length(header2{ii})
              nheaderstr=nheaderstr+1;
          end
          if err
              warning('readbin:read', 'readbin: cell array is truncated: read %d of %d elements', ii-1, nx*ny);
              break;
          end
      end
      if(nheaderstr)
          %if length(header)
          %    header2{end+1}=header;
          %end
          header=header2;
      else
          clear header2;
      end
     case {M_DSP64, M_CSP64}
      nz=zfread(fid,1,'uint64=>double');
      if isempty(nz)
          error('readbin: file is truncated in sparse matrix');
      elseif nz==0 %no index or value arrays follow (see spreaddata in math/spbin.c)
          res=sparse(nx,ny);
          return;
      end
      Jc=zfread(fid,ny+1,'uint64=>uint64');
      Ir=zfread(fid,nz,'uint64=>uint64');
      if magic==M_CSP64
          [P,cnt]=zfread(fid,[2 nz],'double=>double');
      else
          [P,cnt]=zfread(fid,nz,'double=>double');
      end
      convert_sparse=1;
     case {M_DSP32, M_CSP32}
      nz=zfread(fid,1,'uint64=>double');
      if isempty(nz)
          error('readbin: file is truncated in sparse matrix');
      elseif nz==0 %no index or value arrays follow (see spreaddata in math/spbin.c)
          res=sparse(nx,ny);
          return;
      end
      Jc=zfread(fid,ny+1,'uint32=>uint32');
      Ir=zfread(fid,nz,'uint32=>uint32');
      if magic==M_CSP32
          [P,cnt]=zfread(fid,[2 nz],'double=>double');
      else
          [P,cnt]=zfread(fid,nz,'double=>double');
      end
      convert_sparse=1;
     case {M_SSP64, M_ZSP64}
      nz=zfread(fid,1,'uint64=>double');
      if isempty(nz)
          error('readbin: file is truncated in sparse matrix');
      elseif nz==0 %no index or value arrays follow (see spreaddata in math/spbin.c)
          res=sparse(nx,ny);
          return;
      end
      Jc=zfread(fid,ny+1,'uint64=>uint64');
      Ir=zfread(fid,nz,'uint64=>uint64');
      if magic==M_ZSP64 %was M_CSP64: complex single sparse was read as real
          [P,cnt]=zfread(fid,[2 nz],'float=>float');
      else
          [P,cnt]=zfread(fid,nz,'float=>float');
      end
      convert_sparse=1;
     case {M_SSP32, M_ZSP32}
      nz=zfread(fid,1,'uint64=>double');
      if isempty(nz)
          error('readbin: file is truncated in sparse matrix');
      elseif nz==0 %no index or value arrays follow (see spreaddata in math/spbin.c)
          res=sparse(nx,ny);
          return;
      end
      Jc=zfread(fid,ny+1,'uint32=>uint32');
      Ir=zfread(fid,nz,'uint32=>uint32');
      if magic==M_ZSP32 %was M_CSP32: complex single sparse was read as real
          [P,cnt]=zfread(fid,[2 nz],'float=>float');
      else
          [P,cnt]=zfread(fid,nz,'float=>float');
      end
      convert_sparse=1;
     case M_DBL
      [res,cnt]=zfread(fid,[nx,ny],'double=>double');
     case M_FLT
      [res,cnt]=zfread(fid,[nx,ny],'float=>float');
     case M_INT64
      [res,cnt]=zfread(fid,[nx,ny],'int64=>int64');
     case M_INT32
      [res,cnt]=zfread(fid,[nx,ny],'int32=>int32');
     case M_INT16
      [res,cnt]=zfread(fid,[nx,ny],'int16=>int16');
     case M_UINT8
      [res,cnt]=zfread(fid,[nx,ny],'uint8=>uint8');
     case M_CMP
      nread=2*nx*ny;
      [tmp,cnt]=zfread(fid,[2*nx,ny],'double=>double');
      res=squeeze(tmp(1:2:end,:)+1i*tmp(2:2:end,:));
      clear tmp; 
     case M_ZMP
      nread=2*nx*ny;
      [tmp,cnt]=zfread(fid,[2*nx,ny],'float=>float');
      res=squeeze(tmp(1:2:end,:)+1i*tmp(2:2:end,:));
      clear tmp; 
     otherwise
      error('readbin: invalid magic number %d (0x%x)', magic, magic);
    end
    if ~convert_sparse && cnt>=0 && cnt<nread %fread pads [m,n] reads with zeros
        error('readbin: file is truncated: read %d of %d elements', cnt, nread);
    end
    if convert_sparse
        iscmp=any(magic==[M_CSP64 M_CSP32 M_ZSP64 M_ZSP32]);
        if numel(Jc)~=ny+1 || numel(Ir)~=nz || cnt~=nz*(1+iscmp)
            error('readbin: file is truncated in sparse matrix');
        end
        Ic=zeros(nz,1);
        for iy=1:ny
            Ic(Jc(iy)+1:Jc(iy+1))=iy;
        end
        Ir=Ir+1;%start from 1 now.
        P=double(P);
        if iscmp %was size(P,1)==2, which misread real matrices with nz==2
            Pr=P(1,:)';
            Pi=P(2,:)';
            P=Pr+1i*Pi;
        end
        res=sparse(double(Ir),Ic,P,nx,ny);
    end
    