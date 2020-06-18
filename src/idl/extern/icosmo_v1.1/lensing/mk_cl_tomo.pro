;function mk_cl_tomo,cosmo,evol,source,keywords=keywords
;function mk_cl_tomo,fid,cosmo,source
function mk_cl_tomo,fid,cosmo,survey
; Sep 08 - Modified by AnR - check_sv added
; Aug 08 - Modified by AA - header updated
; Aug 08 - Modified by AA suppress underflow error message
; Jul 08 - Modified by AA move source calculation into this routine
; Jul 08 - Modified by AA compatible with new set_fiducial & mk_cosmo
; May 08 - Modified by AA changed to a function
; May 08 - Modified by AA to take in keywords structure
; Mar 08 - Modified by AA fixed bug in z range
; Feb 08 - Modified by AA to speed up
; May 06 - Modified by AA to be double precision
; Feb 2005 - Modified by AR to allow an arbitrary number of redshift bins
; Jan 2005 - Modified by AR to perform tomographics calculation
; July 2003 - Modified by AR to allow the use of different fitting
;             formulae for the non-linear power spectrum
; August 1999 - Written by A. Refregier
;
; ***************// HELP BEGIN //**************
; PURPOSE: compute the lensing power spectra for tomography.
; All auto- and cross-correlation (shear-shear) power spectra for a
; set of redshift bins are computed. In this version, an arbitrary
; number of redshift bins can be treated by turning the source structure
; into an array (see below)
;
; INPUT: fid - fiducial structure created by the routine set_fiducial.pro
;        cosmo - structure containing cosmological
;                information. Created by the routine mk_cosmo.pro
;        survey - survey structure.  Created by mk_survey.pro
;
; OUTPUT: cl structure:
;           n_zbin: number of redshift bins
;           n_cl: number of power spectra (=n_zbin*(n_zbin+1)/2)
;           zbins: indices of the z-bins of each power spectrum
;                  (eg. if zbins(*,2)=[3,4] then cl(*,2) corresponds
;                   to Cl^(34), the cross-correlation PS of bins 3 and 4)
;           l: multipole moment               [1]
;           cl: shear-shear (or kappa-kappa) power spectrum matrix for 
;               each combination of the redshift bins [1]
;----
; Example 1: cl=mk_cl_tomo(fid,cosmo,sv)
; Calculates the lensing power spectra using the inputs set in
; previous routines. e.g. possible call sequence is:
;      fid=set_fiducial('generic')
;      cosmo=mk_cosmo(fid)
;      sv=mk_survey(fid,'sv1')
; ----
; ***************// HELP END //**************

;-- start declarations --
;if not tag_check(keywords,'kl_ran',val=k_ran) then k_ran=[.001d,100000.d]  
;if not tag_check(keywords,'n_k',val=n_k) then n_k=200
;if not tag_check(keywords,'l_ran',val=l_ran) then l_ran=[10.d,2.d4]
;if not tag_check(keywords,'n_l',val=n_l) then n_l=200
;if not tag_check(keywords,'nz',val=n_z) then n_z=50 
;if not tag_check(keywords,'linear',val=linear) then linear=0 
;if not tag_check(keywords,'fit_tk',val=fit_tk) then fit_tk=0b
;if not tag_check(keywords,'fit_nl',val=fit_nl) then fit_nl=2b
 ;-- end declarations --
;k_ran=fid.calc.k_ran
;n_k=fid.calc.n_k
;
;nz=fid.calc.n_z ;!!!!BE CAREFUL ABOUT THIS !!!
;linear=fid.calc.linear
;fit_tk=fid.calc.fit_tk
;fit_nl=fid.calc.fit_nl

;Check if survey is designed for lensing calculations: 
check_sv, survey, 'lens'

current_except=!except
if fid.calc.verbose gt 1 then !except=0 else !except=1
void=check_math()


source=mk_source(fid,cosmo,survey)

l_ran=fid.calc.l_ran
n_l=fid.calc.n_l
k_ran=cosmo.pk.k_ran

l=xgen(l_ran(0),l_ran(1),np=n_l,/log,/double)
z=cosmo.pk.z

n_k=n_elements(cosmo.pk.k)
n_z=n_elements(z)

n_zbin=n_elements(source.z_m)   ; number of z-bins for tomography
n_cl=n_zbin*(n_zbin+1)/2        ; number of power spectra

; form z-bins index vector
zbins=intarr(2,n_cl)            ; z-bin indices for power spectra
k=0
for i=0,n_zbin-1 do begin 
   zbins(*,k)=[i,i]             ; auto-correlations first
   k=k+1
endfor
for i=0,n_zbin-1 do begin       ; then cross-correlations
   for j=i+1,n_zbin-1 do begin
      zbins(*,k)=[i,j]
      k=k+1
   endfor
endfor

if keyword_set(linear) then pk=cosmo.pk.pk_l else pk=cosmo.pk.pk
k_pk=cosmo.pk.k

; compute cosmological distance quantitites
a=interpol(cosmo.evol.a,cosmo.evol.z,z)
sk=interpol(cosmo.evol.sk,cosmo.evol.z,z)
chi=interpol(cosmo.evol.chi,cosmo.evol.z,z)
g=dblarr(n_z,n_zbin)
for i=0,n_zbin-1 do g(*,i)=interpol(source.g(*,i),cosmo.evol.z,z)
;stop
; compute Cl for each l
cl=dblarr(n_l,n_cl)
intgd=dblarr(n_z,n_cl)
k=dblarr(n_z)
jp=where(sk gt 0.d) ; to avoid dividing by zero

for i=0,n_l-1 do begin
   k(*)=-1
   k(jp)=l(i)/(cosmo.const.r0*sk(jp)*cosmo.const.h) ; k [h Mpc^-1]
   intgd=replicate(0.d,n_z,n_cl)
   for j=0,n_z-1 do begin
      if k(j) ge k_ran(0) and k(j) le k_ran(1) then begin
         pk_ij=10.d ^interpol(alog10(pk(*,j)),alog10(k_pk(*)),alog10(k(j)))
         intgd_j=9.d /16.d /cosmo.const.rh^4.d *cosmo.const.omega_m^2.d *cosmo.const.r0*$
                 (1.d /a(j)/sk(j))^2.d *pk_ij(0)/cosmo.const.h^3.d 
         for p=0,n_cl-1 do $
            intgd(j,p)=intgd_j*g(j,zbins(0,p))*g(j,zbins(1,p))
      endif
   endfor
   for p=0,n_cl-1 do cl(i,p)=int_tabulated(chi,intgd(*,p),/double)
endfor

; store in a structure
cl={n_zbin:n_zbin,zbins:zbins,n_cl:n_cl,l:l,cl:cl,l_ran:l_ran}

status = Check_Math() 
if (status and not 32) then message, 'IDL Check_Math() error: ' + StrTrim(status, 2) ; only underflow error are supressed.
!except=current_except
return,cl
end














