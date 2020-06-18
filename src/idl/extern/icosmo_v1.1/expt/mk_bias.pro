function mk_bias,biastype,z,bias_input=bias_input

; Written by Adam Amara 22 July 2008
; PURPOSE: 
; INPUTS:
; OPTIONAL INPUTS: 
; OUTPUTS:
; OPTIONAL OUTPUTS:

n_z=n_elements(z)
bias=make_array(n_z,/double)
if keyword_set(bias_input) then biastype='Input'


case strlowcase(biastype) of 
   'bias1': begin
      bias(*)=sqrt(1.d +z)
   end
   'input': begin
      if not tag_exist(pz_input,'z') then stop  ;!!input pz must contain z
      if not tag_exist(pz_input,'bias') then stop ;!!input pz must contain pz 
      n_zinput=n_elements(pz_input.z)
      bias(*)=interpol(bias_input.bias,bias_input.z,z)
      clip_high=where(z gt max(pz_input.z))
      clip_low=where(z lt min(pz_input.z))
      if (clip_high(0) ne -1) then bias(clip_high)=bias_input(n_zinput-1)
      if (clip_low(0) ne -1) then bias(clip_low)=bias_input(0)     
   end
endcase
return,bias
end


