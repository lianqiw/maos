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

  magic=ulong(1) ;silly int is only 2 byte.
  nx=long64(1)
  ny=long64(1)
  readu,unit,magic,nx,ny
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
        print,nx,ny
        out=ptrarr(nx,ny)
        for iy=0,ny-1 do begin
           for ix=0, nx-1 do begin
              out(ix,iy)=readdata(unit)
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
  return,out
end
