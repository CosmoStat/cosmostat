function func2, z00
; function called by root finding routine fx_root.pro
common pzparams,z0,alpha,beta,zm
numz0=n_elements(z00)
temp=dblarr(numz0)
for i=0,numz0-1 do begin
   z0=z00(i)
   if (z0 le 0.0) then temp(i)=1.0+z0^2 else begin
   norm=qromo('smail',0.d,20.0d,/double)
   temp(i)=(qromo('smail',0.d,zm,/double)/norm) -0.5d
   endelse
endfor
return,temp
end

function smail,z;,z0=z0,alpha=alpha,beta=beta
common pzparams,z0,alpha,beta,zm
pz=z^alpha*exp(-(z/z0)^beta) 
return,pz
end

function mk_smailz0,alphai,betai,zmi
common pzparams,z0,alpha,beta,zm

; Written by Adam Amara 14 August 2008
; PURPOSE: Routine written for calculating z0 of the smail
; distribution when zm, alpha and beta are specified.
; INPUTS: alpha - Smail parameter
;         beta - Smail parameter
;         zm - median redshift
; OUTPUTS: z0 - Smail paramter

if not keyword_set(alphai) then return,'ERROR: must input alpha'
if not keyword_set(betai) then return,'ERROR: must input beta'
if not keyword_set(zmi) then return,'ERROR: must input zm'

alpha=alphai
beta=betai
zm=zmi
z0=newton(zm,'func2',stepmax=0.1,itmax=100.,/double)   

return,z0
end



