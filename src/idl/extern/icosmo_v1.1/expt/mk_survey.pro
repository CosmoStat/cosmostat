function mk_survey,fid,survey,pz_input=pz_input,bias_input=bias_input

; Aug 08 - modified by AA make header compatible with icosmo_help
; Aug 08 - modified by AA general improvements
; Aug 08 - modified by TK to add SNe options
; Aug 08 - modified by AA suppress underflow error message
; Jul 08 - modified by AA - consistent with new set_fiducial
; Modified by AA May 2008 - use keywords structure 
; Written by Adam Amara 23/02/2008
; This routine has been written to replace rd_survey_ideal Written by AR
; ***************// HELP BEGIN //**************
; PURPOSE: produce survey structure with parameters survey
; properties. For the moment the galaxies are distributed using a
; Smail et al distribution.
; 
; INPUTS: fid - fiducial structure set by set_fiducial.pro
;         survey - string specifying the name of the survey of the
;                  experiment. (e.g. 'sv1', 'sv2' and 'sv3'. Make sure
;                  that the chosen experiment contains this survey.)
;
;
; OPTINAL INPUTS: pz_input - user specified probability distribution
;                            of galaxies. Structure of the form
;                            {z:z,pz:pz}. This can either be for all 
;                            galaxies (pz_input.pz is a vector) or
;                            galaxies per slice (pz_input.pz is a 2
;                            dimensional array)
;                 bias_input - User specified galaxy bias. Structure
;                              of the form {z:z,bias:bias}.
;
; OUTPUT: structure containing suvery information
; ----
; Example 1. survey=mk_survey(fiducial,'sv1')
; Calculates the survey structure using the parameters specified in
; fiducial.expt (where the structure fiducial is created by the
; routine set_fiducial.pro).  'sv1' of the experiment is used.
;
; Example 2. survey=mk_survey(fiducial,'sv1',pz_input=pz_input)
; Similar to example 1. However in this example the probablity
; distribution of galaxies (in redshift) is user specified. pz_input
; must be structure containing z and pz (e.g. pz_int={z:z,pz:pz}) 
;
; Example 3. survey=mk_survey(fiducial,'sv2',bias_input=bias_input)
; Similar to example 1. However in this example the galaxy bias is
; user specified. bias must be structure containing z and pz
; (e.g. bias_int={z:z,bias:bias})  
; ----
; ***************// HELP END //**************

current_except=!except
if fid.calc.verbose gt 1 then !except=0 else !except=1
void=check_math()

; *** Check for probe tag: ***
if not keyword_set(survey) then begin
   print,'Error: You have not specified a probe.'
   return,'Error'
endif
void=tag_exist(fid.expt,survey,index=ind)
if not (void) then begin
   print,'The probe you have chosen is not supported.'
   stop
endif

; *** Set probe parameters: ***
a_survey=fid.expt.(ind).a_survey
eff=fid.expt.(ind).eff
sig_int=fid.expt.(ind).sig_int
ng=fid.expt.(ind).ng
z_med=fid.expt.(ind).z_med
zerror=fid.expt.(ind).zerror
n_zbin=fid.expt.(ind).n_zbin
if keyword_set(pz_input) then dndztype='input' else dndztype=fid.expt.(ind).dndztype
dndzp=fid.expt.(ind).dndzp
dndzz=fid.expt.(ind).dndzz
if keyword_set(biasinput) then biastype='input' else biastype=fid.expt.(ind).biastype
ns=fid.expt.(ind).ns
sigmam=fid.expt.(ind).sigmam
delm=fid.expt.(ind).delm
sne_zran=fid.expt.(ind).sne_zran
probes=fid.expt.(ind).probes

; *** make z vector: ***
n_zp=fid.calc.nz_fn
z_ran=fid.calc.ran_z
z=xgen(z_ran(0),z_ran(1),n=n_zp,/double)

; *** Set area and number density of galaxies: ***
a=eff*a_survey/(180.d/!dpi)^2.d ; survey area [srad]
f_sky=a/(4.d*!dpi)              ; fraction of sky covered
n_g=ng*3600.d*(180.d/!dpi)^2.d  ; galaxy surface den. (w/out binning) [srad^-1]

; *** compute the z-distribution specified in set_fiducial: ***
if keyword_set(pz_input) then d
pz=mk_pz(z,dndztype,dndzp,dndzz,n_g=n_g,z_med=z_med,z_min=z_min,z_max=z_max,pz_input=pz_input)
;stop
; *** compute the binned properties: ***
if keyword_set(pz) then begin
   bins= mk_nzbins(pz,z,n_zbin,n_g,zerror=zerror,z_min=z_min,z_max=z_max,z_med=z_med)
endif else begin
   ;rough estimate of the median redshift of the total distribution of galaxies
   if n_elements(z_med) gt 1 then zmed_tot=interpol(z_med,total(n_g,/cum)/total(ng),0.5) else zmed_tot=z_med
   if not keyword_set(z_min) then z_min='NA'
   if not keyword_set(z_max) then z_max='NA'
   ;construct structure:
   bins={z:z,pz_tot:0.d,pz_i:0.d,n_g_i:n_g,z_min:z_min,z_max:z_max,zmean_tot:'NA',zmean_i:'NA',zerror:zerror, zmed_tot:zmed_tot,zmed_i:z_med}
endelse

; *** calculate the bias ***
bias=mk_bias(biastype,z,bias_input=bias_input)

; *** calculate the total number of supernovae ***
nsuper = 0
if ns gt -1 then  nsuper=a_survey*3600.*ns ;ns is number of sne per sq.arcmin

; *** calculate the total number of supernovae redshift bins ***
nbin_ns=round(((delm/sigmam)^2)*1000.) ; assume that the shot noise error is 1000 times smaller than the intrinsic error

; *** create an array of supernovae per bin and redshifts *** 
zbin_sn=0.
nbin_sn=0.
if sigmam gt -1 then begin
	deltaz = (sne_zran[1]-sne_zran[0])/(nbin_ns) ;redshift bin width
	zbin_sn=make_array(nbin_ns)
	nbin_sn=make_array(nbin_ns)
	for i=0,nbin_ns-1 do begin
    		zbin_sn[i]=0
    		nbin_sn[i]=0
    		zbin_sn[i]=sne_zran[0]+(i)*deltaz
    		nbin_sn[i]=nsuper/nbin_ns
	endfor
endif

; *** define the intrinsic velocity dispersion of stars in galaxies for the supernovae *** 
sigmanu=500.0d ;Kms-1	

;*** Construct survey structure: ***
sv={a:a,f_sky:f_sky,sig_int:sig_int,n_g:bins.n_g_i,zmean_bin:bins.zmean_i,zmed_bin:bins.zmed_i,zmean_tot:bins.zmean_tot,zmed_tot:bins.zmed_tot,z:z,pz_tot:bins.pz_tot,pz_bin:bins.pz_i,z_min:bins.z_min,z_max:bins.z_max,survey:survey,dndztype:dndztype,bias:bias,zerror:zerror,sigmam:sigmam,delm:delm,sigmanu:sigmanu,zbin_sn:zbin_sn,nbin_sn:nbin_sn,probes:probes}

status = Check_Math() 
if (status and not 32) then message, 'IDL Check_Math() error: ' + StrTrim(status, 2) ; only underflow error are supressed.
!except=current_except
return,sv
end
