function readbin,file
  openr,unit,file,/get_lun,ERROR = err
  if err ne 0 then begin
     print,'error reading ',file
     return,0
  end
  result=readdata(unit)
  close,unit
  return, result
end
