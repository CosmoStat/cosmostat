function mk_nzbins,pz,z,n_zbin,n_g,zerror=zerror,z_min=z_min,z_max=z_max,z_med=z_med

; Written by Adam Amara 18 July 2008
; PURPOSE: Routine for dividing galaxies into bins 
; INPUTS: pz - input probablity distribution. If this is a vector then
;              this is for all galaxies.  If the pz is 2D array the
;              the routine assumes that the galaxies are already
;              binned.  In this case the total pz is calculated.
;         z - vector of redshifts
;         n_zbin - number of redshift bins
;         n_g - number density of galaxies [#/arcmin^2]. If this is a
;               scalar then this is the total number of galaxies and
;               if this is a scalar (with n_zbin elements) then this
;               is the number density per slice
;         zerror - photoz redshift error
;         z_min - array of the min of the bins
;         z_max - array of the max of the bins
; OUTPUTS: structure containing pz properties

if (n_zbin eq 1) then begin
   zmean=int_tabulated(z,z*pz)/int_tabulated(z,pz)
   zmed=zmed({z:z,p:pz})

   output={z:z,pz_tot:pz,pz_i:pz,n_g_i:n_g,z_min:0,$
           z_max:max(z),zmean_tot:zmean,zmean_i:zmean,zerror:zerror,$
           zmed_tot:zmed,zmed_i:zmed}
   return,output
endif

deltaz=z(1)-z(0)

if not keyword_set(pz) then stop ;!! NEED TO INPUT PZ
dim_pz=(size(pz))(0)
if not ((dim_pz eq 1) or (dim_pz eq 2)) then stop ;!! Dimensions of input pz are wrong
if (n_elements(z) eq n_elements(pz(*,0))) then n_zp = n_elements(z) else stop ;!!pz and z must match

; *** calculating pz of all galaxies: ***
case dim_pz of 
   1: pz_tot=pz
   2:begin
         pz_i_temp=pz
      for i=0,n_zbin-1 do begin
         pz_i_temp(*,i)=pz_i_temp(*,i)/int_tabulated(z,pz_i_temp(*,i),/double)
         pz_i_temp(*,i)=pz_i_temp(*,i)*n_g(i)
      endfor
      pz_tot=total(pz_i_temp,2)
   end
endcase

; *** Calculate zmin and zmax for the slices ***
; conditions for calculating zmin and zmax:
flag1=((not keyword_set(z_min) or not keyword_set(z_max)) and (dim_pz eq 1))
pz_tot=pz_tot/int_tabulated(z,pz_tot,/double)
if flag1 then begin
; compute cumulative distribution function (used for slicing):
   pz_cum=dblarr(n_zp)
   pz_cum=total(pz_tot,/cum)/total(pz_tot)
   z_min=make_array(n_zbin,/double)  ; Minimum boundary of slice
   z_max=make_array(n_zbin,/double)  ; Maximum boundary of slice
   fn_g_i=make_array(n_zbin,/double) ; Fraction of galaxies in slice
   interval=(make_array(n_zbin-1,/index,/double)+1.)/n_zbin
   zbound=interpol(z,pz_cum,interval)
   ; check that the bin widths are not thin compared to z sampling 
   if (n_elements(zbound) gt 2) then begin
      bin_wid=zbound(1:*)-zbound(0:*)
      small_bin=where(bin_wid lt 3.*deltaz)
      if (small_bin(0) eq -1) then begin
         for i=0, n_elements(small_bin)-1 do begin
            zbound(small_bin(i)+1:*)=zbound(small_bin(i)+1:*) +(3.d*deltaz)
         endfor
      endif   
   endif
   z_min=[0.0d,zbound(0:n_zbin-2)]
   z_max=[zbound(0:n_zbin-2),max(z)]
endif

; *** Making slices ***
if (dim_pz eq 2) then begin
   pz_i=pz
endif else begin
   pz_i=make_array(n_zp,n_zbin,/double)
   for i=0,n_zbin -1 do begin
      pz_i(*,i)=mk_cut(pz_tot,z,zerror,z_min(i),z_max(i))
   endfor
endelse

; *** Calculate number of galaxies per slice ***
if (n_elements(n_g) eq 1) then begin
   fn_g_i=make_array(n_zbin,/double)
   for i=0,n_zbin -1 do begin
      fn_g_i(i)=int_tabulated(z,pz_i(*,i),/double)
   endfor
   n_g_i=fn_g_i*n_g/total(fn_g_i)
endif else n_g_i = n_g
; check that the number of elements in n_g_i is correct:
if (n_elements(n_g_i) ne n_zbin) then stop

; *** Calculating the means ***
zmean_tot=int_tabulated(z,z*pz_tot)/int_tabulated(z,pz_tot)
zmean_i=make_array(n_zbin,/double)
for i=0,n_zbin-1 do begin
   zmean_i(i)=int_tabulated(z,z*pz_i(*,i))/int_tabulated(z,pz_i(*,i))
endfor

; *** Calculating medians ***
zmed_tot=zmed({z:z,p:pz_tot})
zmed_i=make_array(n_zbin,/double)
for i=0,n_zbin-1 do begin
   zmed_i(i)=zmed({z:z,p:pz_i(*,i)})
endfor

; *** Construct output structure ***
if not keyword_set(z_min) then z_min = 'Not calculated'
if not keyword_set(z_max) then z_max = 'Not calculated'
if not keyword_set(z_err) then z_err = 'Not used'
output={z:z,pz_tot:pz_tot,pz_i:pz_i,n_g_i:n_g_i,z_min:z_min,z_max:z_max,zmean_tot:zmean_tot,zmean_i:zmean_i,zerror:z_err, zmed_tot:zmed_tot,zmed_i:zmed_i}
;stop
return,output
end

