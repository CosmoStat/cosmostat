function mk_nz,z,dndztype,dndzp,z_m=z_m

; Written by Adam Amara 22 July 2008
; PURPOSE: 
; INPUTS:
; OPTIONAL INPUTS: 
; OUTPUTS:
; OPTIONAL OUTPUTS:

case strlowcase(dndztype) of
   'smail': begin
      alpha=dndzp1
      beta=dndzp2
      if (alpha ne 2) or (beta ne 1.5) then stop ;!! NOT YET POSSIBLE
      z0=z_m/1.412d  ;!!MUST CHANGE!!   ; only valid for alpha=2, beta=1.5
      pz=z^alpha*exp(-(z/z0)^beta) 
   end
   'plane': begin
      ;Construct an array of sv structures:
      z_m_i=dndzp1
      n_g_i=dndzp2
      pz=0b
      if (zerror ne 0.0d) then stop ;The plane option has no z-errors
      if (n_elements(z_m_i) ne n_elements(n_g_i)) then stop;MUST BE SAME SIZE!
   end
   'hist': begin
      zb=dndzp1
      n_g_i=dndzp2
      z_m_i=(zb+zb(1:*))/2.d
      z_min=(shift(zb,1))(1:*)
      z_max=zb(1:*)
   end
   else: begin
      print,'The distribution you have requested in not supported.'
   end
endcase

return,pz

end



