;function mk_cl_cov_tomo,cl,survey,keywords=keywords
function mk_cl_cov_tomo,fid,cl,survey

; Jul 08 - Modified by AA compatible with new set_fiducial & mk_cosmo
; May 08 - Modified by AA changed to a function
; May 08 - Modified by AA to take in keywords structure
; Jul 06 -Modified by AA to include systematics
; May 06 - Modified by AA to be double precision
; Feb 2005 - treatment of arbitrary number of z-bins included by AR
; Jan 2006 - Modified AR to include treatment of tomography with 2 redshift 
;              bins
; Feb 2000 - Written by A. Refregier
;
; PURPOSE: Compute the covariance matrix of the tomographic power spectra
; Both the covraiance per l and errors for bins of l are computed. In this,
; version, an arbitrary number of z-bins can be considered using an array
; of survey structure (se below).
; INPUT: cl: tomographic power spectra from mk_cl_tomo.pro
;        survey: survey parameters structure; set to an array of 
;                structures with each element corresponding a z-bins
; OPTIONAL INPUT: keywords: (as set by the routine keywords.pro). The
;                    following keywords are used:'
;                         l_ran: range of l-value for which to compute the
;                            binned errors (default: [min(l),max(l)]
;                         n_lbin: number of l-bins (default: 20)
; OUTPUT: cl_cov: structure containing:
;             cl_err: 1sigma error for each C_l per l  [1]
;             cl_cov: covariance matrix for the C_ls per l [1]
;             l_b,l_lb,l_hb: central, lower, upper limits of l-bins [1]
;             cl_err_b: binned 1sigma error for C_l [1]
;             cl_cov_b: binned covariance matrix [1]
;             cl_b: C_l at l_b binned values  [1]
; ----
; Example 1: cl_cov=mk_cl_cov_tomo(cl,survey)  
; Computes the covariance matrix using the input cl (calculated using
; mk_cl_tomo) and survey (calculated using mk_survey). Optional
; keywords are set to their inbuilt default values. 
;
; Example 2: cl_cov=mk_cl_cov_tomo(cl,sv,keywords=keywords)  
; Computes the covariance matrix using the input cl (calculated using
; mk_cl_tomo) and survey (calculated using mk_survey). Optional
; keywords are set by the keywords structure (created using set_keywords) 
;
; Example 3: cl_cov=mk_cl_cov_tomo(cl,sv,keywords={l_ran:[100.d,1.d4],n_lbin:40} 
; Computes the covariance matrix using the input cl (calculated using
; mk_cl_tomo) and survey (calculated using mk_survey). The following
; two optional keywords are set: l_ran = [100.d,1.d4],n_lbin= 40
; ----


;-- start declaration --
;if not tag_check(keywords,'l_ran',val=l_ran) then l_ran=[10.d,2.d4]
;if not tag_check(keywords,'n_lbin',val=n_lbin) then n_lbin=20
;-- end declaration --

n_l=n_elements(cl.l)
n_zbin=cl.n_zbin            ; number of redshift bins
n_cl=cl.n_cl    ; number of power spectra
l_ran=cl.l_ran
n_lbin=fid.calc.n_lbin

; check if surveys are compatible
if n_elements(survey.zmed_bin) ne n_zbin then $
  message,'survey array does not have the same number of z-bins as Cl'
if n_elements(uniq(survey.zmed_bin)) lt n_zbin then $
  message,'Warning: z-bins do not appear different in survey structure',/cont

; compute the covariance between each power spectrum
; adapted from Hu & Jain, 2004, PRD 70, 043009,eq. 23
; note: sigma_int is understood to be the rms shear **per component** due to 
;       intrinsic ellipticities and measurement noise
cov_l=dblarr(n_l,n_cl,n_cl)  ; covariance matrix

zbins=cl.zbins
cl_obs=cl.cl              ; observed power spectra (lensing+noise)
for i=0,n_cl-1 do begin
    if zbins(0,i) eq zbins(1,i) then $   ; noise<>0 only for auto-correlations
      cl_obs(*,i)=cl_obs(*,i)+$
                  survey.sig_int^2.d /survey.n_g(zbins(0,i)) 
    ;!!! WARNING: SIG_INT IS Z INDEPENDENT !!!
endfor

geom_l=1.d /(2.d *cl.l+1.d)/survey.f_sky     ; geometric factors


for i=0,n_cl-1 do begin
  for j=i,n_cl-1 do begin
    a1=zbins(0,i) & b1=zbins(1,i)   ; process cov[Cl^(a1,b1),Cl^(a2,b2)]
    a2=zbins(0,j) & b2=zbins(1,j)
    a1a2=where((zbins(0,*) eq a1 and zbins(1,*) eq a2) or $
	       (zbins(1,*) eq a1 and zbins(0,*) eq a2))
    b1b2=where((zbins(0,*) eq b1 and zbins(1,*) eq b2) or $
               (zbins(1,*) eq b1 and zbins(0,*) eq b2))
    a1b2=where((zbins(0,*) eq a1 and zbins(1,*) eq b2) or $
               (zbins(1,*) eq a1 and zbins(0,*) eq b2))
    b1a2=where((zbins(0,*) eq b1 and zbins(1,*) eq a2) or $
               (zbins(1,*) eq b1 and zbins(0,*) eq a2))
    cov_l(*,i,j)=geom_l*(cl_obs(*,a1a2(0))*cl_obs(*,b1b2(0))+$
                         cl_obs(*,a1b2(0))*cl_obs(*,b1a2(0)))
    if i ne j then cov_l(*,j,i)=cov_l(*,i,j)  ; symmetrise the matrix
  endfor
endfor

cov_sys=dblarr(n_l,n_cl,n_cl)  ; covariance matrix
cov_sys_temp=dblarr(n_l,n_cl)

; extract the rms errors from the diagonal
cl_err_l=dblarr(n_l,n_cl)
for i=0,n_cl-1 do cl_err_l(*,i)=sqrt(cov_l(*,i,i))

; compute l-bins for band-averaged errors
l=xgen(l_ran(0),l_ran(1),np=n_lbin,/double)
dlogl=alog(l_ran(1)/l_ran(0))/double(n_lbin)
logl_l=alog(l_ran(0))+findgen(n_lbin)*dlogl
logl_m=logl_l+.5d *dlogl
logl_h=logl_l+dlogl
l_l=exp(logl_l)
l_m=exp(logl_m)
l_h=exp(logl_h)
nl=l_h-l_l

; compute the power spectra at the center of each bin
cl_b=dblarr(n_lbin,n_cl)
for i=0,n_cl-1 do cl_b(*,i)=exp(interpol(alog(cl.cl(*,i)),alog(cl.l),logl_m))

; compute the band average errors for the power spectra
cl_err_b=dblarr(n_lbin,n_cl)
for i=0,n_cl-1 do $
  cl_err_b(*,i)=exp(interpol(alog(cl_err_l(*,i)),alog(cl.l),logl_m))/sqrt(nl)

; Store in structure
cl_cov={cl_err:cl_err_l,cov_l:cov_l,$     ; unbinned quantities
        l_b:l_m,l_lb:l_l,l_hb:l_h,cl_err_b:cl_err_b,cl_b:cl_b}  ; binned quantities

return,cl_cov
end
