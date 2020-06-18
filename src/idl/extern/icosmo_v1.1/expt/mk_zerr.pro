function mk_zerr,pz,z,zerror,zmin,zmax

; Written by Adam Amara 21 July 2008
; PURPOSE: 
; INPUTS:
; OPTIONAL INPUTS: 
; OUTPUTS:
; OPTIONAL OUTPUTS:

pz=pz/int_tabulated(z,pz,/double) ;make sure pz normalised to 1

i1=where(z eq (z(where(z ge zmin)))(0)) 
i2=where(z eq (z(reverse(where(z le zmax))))(0)) 
win=make_array(n_elements(z),/double)
win(i1:i2-1)=1.d
pzout=make_array(n_elements(pz),/double)  
if not keyword_set(zerror) then begin
   pzout=pz*win 
endif else begin
   for j=0,n_elements(z)-1 do begin
      sig_s=zerror*(1+z(j))
      pdf=1.d/(sig_s*sqrt(2.d*!dpi))*exp((-1.d)*(z-z(j))^2.d/2.d/sig_s^2)
      pzout(j)=pz(j)*int_tabulated(z,win*pdf,/double)
   endfor
endelse


return,pzout

end

