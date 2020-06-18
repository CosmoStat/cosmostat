function expt_essence,expt_in

;Written by Adam Amara 18 August 2008
;PURPOSE: This routines contains the survey informations for
;the 'generic' experiment. This is to be used as an example.
;OPTIONAL INPUTS: expt_in - structure containing user specified
;                           inputs.
;OUTPUTS: expt structure that is used in set_fiducial
;
; Experiment probe parameter:
; Lensing:
;    sv_n_zbin - The number of redshift slices (lensing tomography).
;    sv_zerror - z error, sigma(z) = gerror*(1+z).
;    sv_z_med - Median redshift of all galaxies.
;    sv_ng - number density of galaxies [per arc min^2] - this can be
;            a scalr or a vector (where the vector version has the
;            same number of elements at there are tomography bins and
;            the entry in each element is the number of galaxies in
;            that slice)
;    sv_a_survey - area of survey [sq. deg.]
;    sv_eff - masking efficiency
;    sv_sig_int - sigma_e1/G from RRG I
;                   <gamma_1^2>^.5=<gamma_2^2>^.5 from intrinsic+noise    
;    sv_probes - probes that the survey can be used to measure
;    sv_dndzp - parameters of lensing dndz if dndz takes an analytic
;               form (e.g. for Smail et al. p(0) = alpha and (1) = beta)
;    sv_dndzz - parameter 3 of lensing dndz
;    sv_dndztype - form of galaxy distribution for lensing
;    dndztype options: 
;        'smail': galaxies distributed according to the analytic
;                 distribution given by Smail et al
;                  - dndzp - [alpha,beta] 
;                  - dndzz - (optional) The boundary of the bins. This
;                    should be a vector where the number of elements
;                    is sv_n_zbin+1.  eg for 3 bins we could have
;                    dndzz=[0.d,0.5d,1.0d,3.0]. If dndzz is not set
;                    (i.e. dndzz=0) then the bins will be made so that
;                    each bin contains (roughly) the same number of
;                    galaxies. 
;                  - ng - should be scalar (i.e. total galaxies in the
;                    survey). The number of galaxies per bin will be
;                    calculated by mk_survey
;                               
;        'plane': galaxies are placed on a planes (or sheets) in redshift
;                  - dndzp - not used
;                  - dndzz - redshift of the (scalar or vector) 
;                  - ng - number of galaxies [# gals/arcmin^2] per
;                    plane.  number of elements should match those of dndzz 
;                  - z_m, z_error and n_zbins should not be set
;
;        'hist': galaxies are distributed in histrograms
;                  - dndzp - not used
;                  - dndzz - The boundary of the bins. This
;                    should be a vector where the number of elements
;                    is sv_n_zbin+1.  eg for 3 bins we could have
;                    dndzz=[0.d,0.5d,1.0d,3.0]
;                  - ng - number of galaxies [# gals/arcmin^2] per
;                    plane.  The number of elements should be = number
;                    of bins
;
;  Not yet implimented:'Gauss'
; 
;  References: 0510026 Sollerman et al. (2005); 0701043 Miknaitis et
;  al. (2007)
                    
names=['wide','spectra','deep']  
  
; Deep Survey:
sv3_n_zbin=10                   ;# of z bins
sv3_zerror=0.04d                ;z error 
sv3_z_med=1.5d                  ;median z
sv3_ng=100.0d                   ;# galaxies [# gals/arcmin^2]
sv3_a_survey=2.d                ;area [sq. deg.]
sv3_eff=1.0d                    ;masking eff.
sv3_sig_int=-1                  ;intrinsic+noise
sv3_dndztype='smail'            ;type of gal dist
sv3_dndzp=[2.0d,1.5d]           ;param of gal dist
sv3_dndzz=0                     ;z info for gal dist
sv3_biastype='bias1'            ;bias type
sv3_ns=2.8d-2                   ;# Sne [# SNe/arcmin^2]
sv3_sigmam=0.15d                ;intrisic scatter in Sne mag [magnitude] 
sv3_delm=0.05d                  ;magnitude precision [magnitude]
sv3_sne_zran=[0.2d,0.8d]        ;Sne redshift range
sv3_probes=['sne']              ;possible probes 

;*** User Customized (sv3): ***
if tag_check(expt_in,'sv3_n_zbin',val=val) then sv3_n_zbin=val 
if tag_check(expt_in,'sv3_zerror',val=val) then sv3_zerror=val      
if tag_check(expt_in,'sv3_z_med',val=val) then sv3_z_m=val
if tag_check(expt_in,'sv3_ng',val=val) then sv3_ng=val
if tag_check(expt_in,'sv3_a_survey',val=val) then sv3_a_survey=val
if tag_check(expt_in,'sv3_eff',val=val) then sv3_eff=val
if tag_check(expt_in,'sv3_sig_int',val=val) then sv3_sig_int=val
if tag_check(expt_in,'sv3_dndztype',val=val) then sv3_dndztype=val
if tag_check(expt_in,'sv3_dndzp',val=val) then sv3_dndzp=val
if tag_check(expt_in,'sv3_dndzz',val=val) then sv3_dndzz=val
if tag_check(expt_in,'sv3_probes',val=val) then sv3_probes=val
if tag_check(expt_in,'sv3_ns',val=val) then sv3_ns=val
if tag_check(expt_in,'sv3_sigmam',val=val) then sv3_sigmam=val
if tag_check(expt_in,'sv3_delm',val=val) then sv3_delm=val
if tag_check(expt_in,'sv3_sne_zran',val=val) then sv3_sne_zran=val

sv3= {n_zbin:sv3_n_zbin,zerror:double(sv3_zerror),z_med:double(sv3_z_med),ng:double(sv3_ng),a_survey:double(sv3_a_survey),eff:double(sv3_eff),sig_int:double(sv3_sig_int),dndztype:sv3_dndztype,dndzp:double(sv3_dndzp),dndzz:double(sv3_dndzz),biastype:sv3_biastype,ns:double(sv3_ns),sigmam:double(sv3_sigmam),delm:double(sv3_delm),sne_zran:double(sv3_sne_zran),name:names(2),probes:sv3_probes}

;expt=create_struct('surveys',names,names(2),sv3)
expt={sv3:sv3}

return,expt

end

