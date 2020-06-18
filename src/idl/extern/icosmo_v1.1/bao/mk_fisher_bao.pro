;function mk_fisher_bao, cl_info=cl_info;, sv, cl_info=cl_info
function mk_fisher_bao,fid,sv,bao_info=bao_info,cosmo=cosmo

;;=============================================================
;;TEST MODE
;;To be replaced by a call to mk_survey
;cmrestore,filename= 'example_survey.sav'
;;mk_survey should also outputt bias_bin
;;linterp,sv.z, sv.bias, (sv.z_min+sv.z_max), bias_bin
;bias_bin=interpol(sv.bias,sv.z,(sv.z_min+sv.z_max)/2.d)
;sv = create_struct(sv, 'bias_bin', bias_bin)
;;END TEST MODE
;;=============================================================

; Aug 08 - modified by AnR to output eigenvalues and do matrix checks in double precision
; Jul 08 - modified by AA make compatible with v0.11
; Written by Anais Rassat, July 2008
; July 2008, Updated by Anais Rassat. Calculates derivatives with
;            4-point function.
; Structure similiar to mk_fisher_lens written by Adam Amara July 2008.
; PURPOSE: This routine calculates the bao Fisher matrix using the
; Formalism of Blake et al. 2006 and Parkinson et al. 2006.
;
; INPUTS: fiducial - Structure with the cosmological model
; OPTIONAL INPUTS:

; -- Start: check sv can do lensing --
void=where(sv.probes eq 'bao',count)
if (count eq 0) then begin
   print,'Error:This survey cannot be used for the current implimentation of BAOs'
   return,-1
endif

;if (sv.dndztype ne 'hist') then begin
;   print,'Warning: The BAO calcualtion is strictly only valid for hist mode and needs more testing.'
;   return,-1
;endif

; -- End: check sv can do lensing --

; -- Start: Define central model --
;  fiducial= set_fiducial({ran_z:[1e-3,3.d0]})
  cosmo = mk_cosmo(fid)      ; central cosmology
  bao = mk_bao(cosmo,sv)
  bao_cov=mk_bao_cov(cosmo,bao,sv)
; -- End: Define central model --

  
; -- Start: Define Fisher matrix parameters --
; construct vector of central values parameters:
  p=[cosmo.const.omega_m,cosmo.const.w0,cosmo.const.wa,cosmo.const.h,cosmo.const.omega_b,$
     cosmo.const.omega_l]
  pname=['!7X!6!im!n','w0','wa','h','!7X!6!ib!n','!7X!6!il!n']
;Check that the number of parameters of the two parameter vectors match
  if (n_elements(p) ne n_elements(pname)) then begin
     print,'The number of parameters in pname and p vectors do not match.'
     return,'Error: # of parameters.'
  endif else n_p=n_elements(p)

;print cosmology parameters
; -- End: Define Fisher matrix parameters --

; -- Start: Calculate derivatives --
  n_z = n_elements(sv.z_min)
; array that will contain dy/dp and dyp/dp functions
  dydp  = dblarr(n_p, n_z)      ;d(r(z)/s)/dp
  dypdp = dblarr(n_p, n_z)      ;d(r'(z)/s)/dp
  
; vary omega_m
  if keyword_set(fid.calc.speed) then $
     dbaodp_temp=mk_dbaodp_2p(fid,sv,'omega_m') else $
        dbaodp_temp=mk_dbaodp_4p(fid,sv,'omega_m')
  
  dydp(0,*)=dbaodp_temp.dydp_bin
  dypdp(0,*)=dbaodp_temp.dyprimedp_bin

; vary w_n
  if keyword_set(fid.calc.speed) then $
     dbaodp_temp=mk_dbaodp_2p(fid,sv,'w0') else $
        dbaodp_temp=mk_dbaodp_4p(fid,sv,'w0')
  
  dydp(1,*)=dbaodp_temp.dydp_bin
  dypdp(1,*)=dbaodp_temp.dyprimedp_bin

; vary w_a
  if keyword_set(fid.calc.speed) then $
     dbaodp_temp=mk_dbaodp_2p(fid,sv,'wa') else $
        dbaodp_temp=mk_dbaodp_4p(fid,sv,'wa')
  
  dydp(2,*)=dbaodp_temp.dydp_bin
  dypdp(2,*)=dbaodp_temp.dyprimedp_bin
  
; vary h
  if keyword_set(fid.calc.speed) then $
     dbaodp_temp=mk_dbaodp_2p(fid,sv,'h') else $
        dbaodp_temp=mk_dbaodp_4p(fid,sv,'h')
  
  dydp(3,*)=dbaodp_temp.dydp_bin
  dypdp(3,*)=dbaodp_temp.dyprimedp_bin

; vary sigma_8
;  dydp_temp=mk_dydp_4p(fiducial,cosmo,'sigma8', survey=sv)
;  dypdp_temp=mk_dyprimedp_4p(fiducial,cosmo,'sigma8', survey=sv)
;  
;  dydp(4,*)=dydp_temp
;  dypdp(4,*)=dypdp_temp

; vary omega_b
  if keyword_set(fid.calc.speed) then $
     dbaodp_temp=mk_dbaodp_2p(fid,sv,'omega_b') else $
        dbaodp_temp=mk_dbaodp_4p(fid,sv,'omega_b')
  
  dydp(4,*)=dbaodp_temp.dydp_bin
  dypdp(4,*)=dbaodp_temp.dyprimedp_bin

; vary n
;  dydp_temp=mk_dydp_4p(fiducial,cosmo,'n', survey=sv)
;  dypdp_temp=mk_dyprimedp_4p(fiducial,cosmo,'n', survey=sv)
;  
;  dydp(6,*)=dydp_temp
;  dypdp(6,*)=dypdp_temp

if keyword_set(fid.cosmo.curv) then begin
; vary omega_l
   if keyword_set(fid.calc.speed) then $
      dbaodp_temp=mk_dbaodp_2p(fid,sv,'omega_l') else $
         dbaodp_temp=mk_dbaodp_4p(fid,sv,'omega_l')
   
   dydp(5,*)=dbaodp_temp.dydp_bin
   dypdp(5,*)=dbaodp_temp.dyprimedp_bin

endif
;-- End: Calculate derivatives --


; For Flat cosmology need to remove Omega_v entries
if not keyword_set(fid.cosmo.curv) then begin
   keep = make_array(n_p,value=1)
   ind=where(pname eq '!7X!6!il!n')
   keep(ind) = -666
   ii=where(keep eq 1)
   pname=pname(ii)
   dydp=dydp(*,ii,*)
   dypdp=dypdp(*,ii,*)
   n_p=n_p-1    
endif

; compute fisher matrix
if (fid.calc.verbose ge 1) then print,'Computing Fisher Matrices'
f=dblarr(n_p,n_p)          ; Fisher matrix for n_p parameters
f[*,*]=0.d0
for i=0,n_p-1 do begin
   for j=0,n_p-1 do begin
      f(i,j) = total(double(dypdp(i,*)*dypdp(j,*)/bao_cov.delta_rad^2))+total(double(dydp(i,*)*dydp(j,*)/bao_cov.delta_trans^2)) 
   endfor
endfor


;;Checking if the diag. elements of Fisher matrix are positive
;check_diag = check_diag(f)
;if check_diag eq 1 then stop
;;Checking if the eigenvalues are positive
;check_eig  = check_eig(f, d=evalues, eps=eps)       ; can use keyword double he;re
;if check_eig eq 2 then stop

;Creating Fisher structure
cov=la_invert(f,/double)
sigma_fix=dblarr(n_p) & sigma_marg=dblarr(n_p)
sigma_fix_diag=dblarr(n_p) & sigma_marg_diag=dblarr(n_p)
sigma_fix_auto=dblarr(n_p) & sigma_marg_auto=dblarr(n_p)
for i=0,n_p-1 do begin
   sigma_fix(i)=sqrt(1./f(i,i)) & sigma_marg(i)=sqrt(cov(i,i))
endfor


if keyword_set(printit) then begin
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

;construct a structure with extra bao information:
bao_info={cosmo:cosmo,bao:bao,dydp:dydp, dypdp:dypdp, pname:pname,bao_cov:bao_cov}

;construct the fisher matrix:
fisher={pname:pname(0:n_p-1),p:p(0:n_p-1),f:f,cov:cov,sigma_marg:sigma_marg,sigma_fix:sigma_fix,status:status, evalues:evalues}

return,fisher

end
