function mk_pz,z,dndztype,dndzp,dndzz,n_g=n_g,z_med=z_med,z_min=z_min,z_max=z_max,pz_input=pz_input
  
; Modified by Anais Rassat August 2008. Fixed slight bug relating to dndzz for
; smail case
; Written by Adam Amara 22 July 2008
; PURPOSE: Routine for calculating the properties of the total galaxy
; distribution. 
; INPUTS: z - vector for redshifts
;         dndztype - type of distribution: Smail, plane, hist and
;                    input
;         dndzp - parameters of the distribution (e.g. for Smail et
;                 al these are alpha and beta)
; OPTIONAL INPUTS: n_g - # of galaxies [#/arcmin^2]
;                  z_med - median redshift
;                         dimensions as zmin)
;                  pz_input - structure containing manually input z
;                             and pz. These will be resampled to match
;                             the z used in mk_cosmo 
; OUTPUTS: pz - structure containing 
; OPTIONAL OUTPUTS:zmin - vector with minimum boundary for z cut
;                  zmax - vector with max boundarf of z cut (same


if keyword_set(pz_input) then dndztype='Input'

case strlowcase(dndztype) of
   'smail': begin
      alpha=dndzp(0)
      beta=dndzp(1)
      if not keyword_set(z_med) then begin
         print,'median redshift must be specified for Smail distibution'
         stop
      endif
      if (max(z) lt 3.*z_med) then begin
         print,'Warning: the redshift range you have chosen is not large enough for this p(z) distibution'
         stop
      endif
      z0=mk_smailz0(alpha,beta,z_med) 
      pz=z^alpha*exp(-(z/z0)^beta) 
      if keyword_set(dndzz) then begin
         if (dndzz[0] ne -1) then begin
            zb=dndzz
            z_min=(shift(zb,1))(1:*)
            z_max=zb(1:*)     
         endif
      endif
   end
   'plane': begin
      if (n_elements(dndzz) ne n_elements(n_g)) then stop ;MUST BE SAME SIZE!
      if keyword_set(zerror) then begin
         if (zerror ne 0.0d) then stop ;The plane option has no z-errors
      endif
      z_med=dndzp
      pz=0b
   end
   'hist': begin
      zb=dndzz
      if (n_elements(zb) ne n_elements(n_g)+1) then stop ;sizes must match
      n_zbin=n_elements(zb)-1
;      z_m_i=(zb+zb(1:*))/2.d
      z_min=(shift(zb,1))(1:*)
      z_max=zb(1:*)
      z_med=(z_max+z_min)/2.0d
      pz=0b
;      pz=make_array(n_elements(z),n_zbin,/double,val=0.0d)
;      for i=0,n_zbin-1 do begin
;         j1=where(z eq (z(where(z ge z_min(i))))(0)) 
;         j2=where(z eq (z(reverse(where(z le z_max(i)))))(0)) 
;         pz(j1:j2-1,i)=1.d /(z_max(i) - z_min(i))
;      endfor
   end
   'input': begin
      if not tag_exist(pz_input,'z') then stop  ;!!input pz must contain z
      if not tag_exist(pz_input,'pz') then stop ;!!input pz must contain pz 
      size_pz=size(pz_input.pz)
      if (size_pz(0) eq 2) then n_pz=size_pz(2) else n_pz=1
      pz=make_array(n_elements(z),n_pz,/double,val=0.0d)
      for i=0,n_pz-1 do begin
         pz_temp=interpol(pz_input.pz,pz_input.z,z) 
         clip_high=where(z gt max(pz_input.z))
         clip_low=where(z lt min(pz_input.z))
         if (clip_high(0) ne -1) then pz_temp(clip_high)=0.0d
         if (clip_low(0) ne -1) then pz_temp(clip_low)=0.0d
         pz(*,i)=pz_temp
      endfor
      if keyword_set(dndzz) then begin
         zb=dndzz
         z_min=(shift(zb,1))(1:*)
         z_max=zb(1:*)     
      endif
   end
   else: begin
      print,'The distribution you have requested in not supported.'
   end
endcase

return,pz

end



