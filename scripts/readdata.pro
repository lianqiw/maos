;; 2010-11-08: Initial version
;; 2010-12-13: Contributed by Mark Chun on 2010-11-29 to make it work for cell arrays and headers

function readdata,unit

  M_CSP64=25600
  M_SP64=25601
  M_CSP32=25606
  M_SP32=25607
  M_DBL=25602
  M_INT64=25603
  M_CMP=25604
  M_INT32=25605
  MC_CSP=25616
  MC_SP=25617
  MC_DBL=25618
  MC_INT64=25619
  MC_CMP=25620
  MC_INT32=25621
  MCC_ANY=25633
  MCC_DBL=25634
  MCC_CMP=25636
  MAT_SP=65281
  MAT_CSP=65282
  M_HEADER=25856    ;; MC:  2010-11-29 Added M_HEADER code to handle header text

  magic=ulong(1) ;silly int is only 2 byte. we need 4 byte int
  nx=long64(1) ;we need 8 byte int.
  ny=long64(1)
  readu,unit,magic
  header='' ;; MC:  2010-11-29 Added M_HEADER code to handle header text
  while (magic EQ M_HEADER) DO BEGIN
        nlen = ulong64(1)
        readu, unit, nlen
        achar = byte(1)
        FOR i=0,nlen-1 DO BEGIN
           readu, unit, achar
           header = header+STRING(achar)
        ENDFOR
        PRINT, header
        nlen2 = ulong64(1)
        readu, unit, nlen2
        magic2=1UL
        readu, unit, magic2
        IF ( nlen NE nlen2 OR magic NE magic2 ) THEN STOP
        readu, unit, magic
  endwhile           ;; MC:  2010-11-29 Added M_HEADER code to handle header text
  readu, unit, nx,ny
  print, magic, nx, ny
  IF ( nx NE 0 and ny NE 0 ) THEN BEGIN
   switch magic of
     MCC_ANY:
     MCC_DBL:
     MCC_CMP:
     MC_CSP:
     MC_SP:
     MC_DBL:
     MC_CMP:
     MC_INT32:
     MC_INT64: begin
        out=ptrarr(nx,ny,/ALLOCATE_HEAP)  ;; MC: Added /ALLOCATE_HEAP
        for iy=0,ny-1 do begin
           for ix=0, nx-1 do begin
              *out(ix,iy)=readdata(unit)  ;; MC: Changed this to *(out)
           end
        end
        break
     end
     M_SP64:
     M_SP32:
     M_CSP64:
     M_CSP32:begin
        break
     end
     M_DBL:begin
        out=dblarr(nx,ny)
        readu,unit,out
        break
     end
     M_INT64:begin
        out=long64arr(nx,ny)
        readu,unit,out
        break
     end
     M_INT32:begin
        out=intarr(nx,ny)
        readu,unit,out
        break
     end
     M_CMP:begin
        out=dcomplexarr(nx,ny)
        readu,unit,out
        break
     end
   end
  ENDIF else out = -1   ;; MC: added setting out=-1 if no match
  
  return,out
end
