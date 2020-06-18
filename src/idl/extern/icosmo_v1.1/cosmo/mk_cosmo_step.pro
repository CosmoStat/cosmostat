function mk_cosmo_step,fiducial,step_param,step_size

; Written by Adam Amara 15th May 2008
; PURPOSE: This routine can be used to make a small step in one of the
; cosmological parameters. This is useful when doing a Fisher matrix
; calculation.
; INPUT: fiducial: structure containing the fiducial cosmology
;                  parameters
;        step_param: name of parameter to be changes (string)
;        step_size: size of the step to be taken
; ----

cosmo_temp=fiducial

if not keyword_set(step_size) then step_size=0.006d

if not tag_exist(cosmo_temp,step_param) then begin
   print,'Error: The parameter you are trying to change is not one of the main model parameters.'
   return,'Error: not a parameter'
endif

case step_param of
   'omega_m': begin
      if (cosmo_temp.curv eq 0) then begin
         cosmo_temp.omega_m= cosmo_temp.omega_m+step_size
         cosmo_temp.omega_l= 1.d - cosmo_temp.omega_m
      endif else cosmo_temp.omega_m= cosmo_temp.omega_m+step_size
   end
   'omega_l': begin
      if (cosmo_temp.curv eq 0) then begin
         cosmo_temp.omega_l= cosmo_temp.omega_l+step_size
         cosmo_temp.omega_m= 1.d - cosmo_temp.omega_l
      endif else cosmo_temp.omega_l= cosmo_temp.omega_l+step_size
   end
   else: begin
      void=tag_exist(cosmo_temp,step_param,index=ind)
      cosmo_temp.(ind)= cosmo_temp.(ind) +step_size
   end
endcase
return,cosmo_temp

end




