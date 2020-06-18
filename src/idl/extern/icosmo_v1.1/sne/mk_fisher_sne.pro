function mk_fisher_sne,fid,sv,cosmo=cosmo

; Written by Tom Kitching, August 2008
; Structure similiar to mk_fisher_lens written by Adam Amara July 2008.
; PURPOSE: This routine calculates the supernova type-IA 
; Fisher matrix using the
; Formalism of Tegmark et al. (1998) and Huterer & Turner (2001).
;
; INPUTS: fiducial - Structure with the cosmological model
; OPTIONAL INPUTS:

; -- Start: check sv can do sne --
void=where(sv.probes eq 'sne',count)
if (count eq 0) then begin
   print,'Error:This survey cannot be used for the current implimentation of Supernovae'
   return,-1
endif

; -- End: check sv can do sne --

; -- Start: Define central model --
if not keyword_set(cosmo) then  cosmo=mk_cosmo(fid,/nopk)
  sne_cov= mk_sne_cov(cosmo,sv) 
; -- End: Define central model --
  
; -- Start: Define Fisher matrix parameters --
; construct vector of central values parameters:
  p=[cosmo.const.omega_m,cosmo.const.w0,cosmo.const.wa,cosmo.const.h,$
	cosmo.const.omega_l]
  pname=['!7X!6!im!n','w0','wa','h','!7X!6!il!n']
; Check that the number of parameters of the two parameter vectors match
  if (n_elements(p) ne n_elements(pname)) then begin
     print,'The number of parameters in pname and p vectors do not match.'
     return,'Error: # of parameters.'
  endif else n_p=n_elements(p)

;print cosmology parameters
; -- End: Define Fisher matrix parameters --

; -- Start: Calculate derivatives --

nz_int=n_elements(sv.zbin_sn)

; array that will contain dCl/dp functions
ddldp=dblarr(nz_int,n_p) 

; vary omega_m
ddldp_temp=mk_ddldp_4p(fid,sv,'omega_m')
ddldp(*,0) =ddldp_temp

; vary w_n
ddldp_temp=mk_ddldp_4p(fid,sv,'w0')
ddldp(*,1) =ddldp_temp

; vary w_a
ddldp_temp=mk_ddldp_4p(fid,sv,'wa')
ddldp(*,2) =ddldp_temp

; vary h
ddldp_temp=mk_ddldp_4p(fid,sv,'h')
ddldp(*,3) =ddldp_temp

if fid.cosmo.curv eq 1 then begin
; vary omega_l
    ddldp_temp=mk_ddldp_4p(fid,sv,'omega_l')
    ddldp(*,4) =ddldp_temp
endif
;-- End: Calculate derivatives --

; For Flat cosmology need to remove Omega_v entries
if fid.cosmo.curv ne 1 then begin
   keep = make_array(n_p,value=1)
   ind=where(pname eq '!7X!6!il!n')
   keep(ind) = -666
   ii=where(keep eq 1)
   pname=pname(ii)
   ddldp=ddldp(*,ii,*)
   n_p=n_p-1    
endif

; compute fisher matrix
if fid.calc.verbose ge 1 then begin
    print,'Computing Fisher Matrices'
endif
f=dblarr(n_p,n_p)          ; all power spectra with full covariance
for i=0,n_p-1 do begin
  for j=0,n_p-1 do begin
      intgd=0.
      for k=1,nz_int-1 do begin
          intgd=intgd+transpose(reform(ddldp(k,i)))#$
            invert(reform(sne_cov.cov_z(k)))#reform(ddldp(k,j))
      endfor
      f(i,j)=f(i,j)+intgd    ;int_tabulated(evol.z,intgd,/double)
  endfor
endfor

cov=la_invert(f,/double)
sigma_fix=dblarr(n_p) & sigma_marg=dblarr(n_p)
for i=0,n_p-1 do begin
  sigma_fix(i)=sqrt(1./f(i,i)) & sigma_marg(i)=sqrt(cov(i,i))
endfor

if fid.calc.verbose ge 1 then begin
  print,'Fisher matrix: F_ij:',f
  print,'sigma (fixed):',sigma_fix
  print,'cov(pi,pj):',cov
  print,'sigma (marginalised):',sigma_marg
endif

;check that the fisher and covariance matrices have positive diagonals
;and eigen values:
f_check=check_matrix(f, d=f_ev, /double)
cov_check=check_matrix(cov, d=c_ev,/double)
status={f_d:f_check.CHECK_D,f_ev:f_check.CHECK_ev,$
        cov_d:cov_check.CHECK_D,cov_ev:cov_check.CHECK_ev}
evalues={fish:f_ev, cov:c_ev}

fisher={pname:pname(0:n_p-1),p:p(0:n_p-1),f:f,cov:cov,sigma_marg:sigma_marg,sigma_fix:sigma_fix}

return,fisher

end
