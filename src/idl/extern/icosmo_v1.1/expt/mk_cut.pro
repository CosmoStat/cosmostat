function mk_cut,pz,z,zerror,zmin,zmax

; Aug 08 - modified by AA
; Written by Adam Amara 21 July 2008
; PURPOSE: This routine takes a total galaxy distribution and makes
; cuts in photoz space at zmin and zmax.  This produces a pz for in
; spectra space that includes photoz errors.
; INPUTS: pz - distribution of all galaxies (vector)
;         z - redshift (vector)
;         zerror - photoz error ; sig(z) = zerror*(1+z)
;         zmin - zmin at which slicing is done
;         zmax - zmax at which slicing is done
; OUTPUTS: returns pz after slicing
;

pz=pz/int_tabulated(z,pz,/double) ;make sure pz normalised to 1
deltaz=z(1)-z(0)
if (zerror lt deltaz/2.) then zerr=0 else zerr=zerror

pzout=make_array(n_elements(pz),/double)  

if not keyword_set(zerr) then begin
   ind=where((z gt zmin) and (z lt zmax))
   win=make_array(n_elements(z),/double,val=0.0d)
   win(ind)=1.d
   pzout=pz*win 
endif else begin
jmin=max([(zmin-(zerr*(1+zmin)*10.))/deltaz,0.d])
jmax=min([(zmax+(zerr*(1+zmax)*10.))/deltaz,n_elements(z)-1])
ztemp=[xgen(zmin,zmax,np=100.,/double),zmax]
pdftemp=make_array(n_elements(ztemp),/double,val=0.0d)
   for j=jmin,jmax do begin
      pdftemp(*)=0.0d
      sig_s=zerr*(1+z(j))
      pdftemp=1.d/(sig_s*sqrt(2.d*!dpi))*exp((-1.d)*(ztemp-z(j))^2.d/2.d/sig_s^2)
      pzout(j)=pz(j)*int_tabulated(ztemp,pdftemp,/double)
   endfor
endelse
return,pzout

end

