function mk_source,fid,cosmo,survey
; Modified by AA Aug 2008 - making function compatible with v0.10.4
; Written by Adam Amara 24th Feb 2008
; This is a modified version of mk_source writen by AR (most of which
; has been moved to mk_srouce_one)
;
; PURPOSE: Calculates an array of source structures. The length of the
; array is the number of redshift slices, which is set in
; rd_survey. Each individual element of the array is set by the
; routine mk_source_one. See this routine for more details.
; INPUT: cosmo: cosmology structure read by rd_cosmo.pro
;        evol: evol structure made by mk_evol.pro
; OPTIONAL INPUT: survey: survey parameter structure from
;                         mk_survey.pro. These will be superseded by
;                         the keywords. 
;                 !!! WARNING: KEYWORDS ARE NOT BEING USED !!
;                 !!! WARNING: MUST USE SURVEY STRUCTURE  !!
;                 keywords: keywords structure (see mk_keywords for
;                           more details). The following keywords are
;                           used by  mk_source and mk_source_one:
;                 !!! WARNING: KEYWORDS ARE NOT BEING USED !!
;                 !!! WARNING: MUST USE SURVEY STRUCTURE  !!
;
; OUTPUT: Returns an array containing source structures. (all arrays
; indexed by evol tables) 
;            pz: Pz(z): probablility of finding a source at resdhift z
;                       normalised so that int dz Pz =1
;            pchi: Pchi(chi): probablility of finding a source at comoving
;                       radius chi, normalised so that int dchi Pchi =1 
;                       (Pchi dchi=Pz dz)
;            g: lensing radial weight function for Smail et al. distribution
;                        in units of R_0 [1]
;            g1: lensing radial weight function for sources on a sheet
;                       at zs=z_median (z0: parameter for Smail et al.
;                       distribution) in units of R_0 [1]
; ----
; Example 1: source=mk_source(cosmo,evol,survey=sv)
; Calculates source structure using the parameters set in the survey
; structure. 
;
; Example 2: source=mk_source(cosmo,evol)
; Calculates the source structure using the built in default values -
; z_m=1.d, alpha=2.d, beta=1.5d, analytic=1 and one redshift slice is
; used. 
;
; Example 3: source=mk_source(cosmo,evol,keywords=keywords)
; Calculates the source structure using the values set in the keywords
; structure. !!! WARNING: NOT YET IMPLEMENTED !!
; ----

;if keyword_set(survey) then begin
;   for i=0,n_elements(survey)-1 do begin
;      mk_source_one,cosmo,evol,source_i,survey=survey(i),keywords=keywords 
;      if (i eq 0) then source=source_i else source=[source,source_i]
;   endfor
;endif else mk_source_one,cosmo,evol,source,keywords=keywords

n_z=n_elements(cosmo.evol.z)
if not keyword_set(cosmo) then stop ;!!need this structure
if not keyword_set(survey) then stop ; !! need this structure

z_m=survey.zmed_bin             ; compute the median redshift
n_zbins=n_elements(z_m)
pz=make_array(n_z,n_zbins,/double)
pchi=make_array(n_z,n_zbins,/double)
for i=0,n_zbins-1 do begin
   pz(*,i)=interpol(survey.pz_bin(*,i),survey.z,cosmo.evol.z) ; interpolate on output grid ;!!BETTER WAY OF DOING THIS??
endfor
pz=pz>0.d ; make sure pz is +ve

; convert p_z to p_chi
chi_s=interpol(cosmo.evol.chi,cosmo.evol.z,z_m)
dzdchi= 1.d /cosmo.const.rh*cosmo.evol.hc/cosmo.const.h0 ; valid in general
for i=0,n_zbins-1 do begin
   pchi(*,i)=pz(*,i)*dzdchi
endfor

; normalize so that \int dchi p_chi =1
for i=0,n_zbins-1 do begin
   norma=int_tabulated(cosmo.evol.chi,pchi(*,i),/double)*cosmo.const.r0
   pchi(*,i)=pchi(*,i)/norma
   pz(*,i)=pz(*,i)/norma
endfor
g=dblarr(n_z,n_zbins)
g1=dblarr(n_z,n_zbins)

case 1 of
   cosmo.const.omega eq 1.: begin
      for j=0,n_zbins-1 do begin
         g1(*,j)=2.*cosmo.evol.chi*(chi_s(j)-cosmo.evol.chi)/chi_s(j)
         for i=1,n_z-2 do begin
            g(i,j)=2.*cosmo.const.r0* $
                   int_tabulated(cosmo.evol.chi(i:n_z-1),pchi(i:n_z-1,j)*$
                                 cosmo.evol.chi(i)*$
                                 (cosmo.evol.chi(i:n_z-1)-cosmo.evol.chi(i))/$
                                 cosmo.evol.chi(i:n_z-1),/double)
         endfor
      endfor
   end
   cosmo.const.omega gt 1.: begin
      for j=0,n_zbins-1 do begin
         g1(*,j)=2.*sin(cosmo.evol.chi)*sin(chi_s(j)-cosmo.evol.chi)/sin(chi_s(j))
         for i=1,n_z-2 do begin
            g(i,j)=2.*cosmo.const.r0*$
                   int_tabulated(cosmo.evol.chi(i:n_z-1),pchi(i:n_z-1,j)*$
                                 sin(cosmo.evol.chi(i))*$
                                 sin(cosmo.evol.chi(i:n_z-1)-cosmo.evol.chi(i))/$
                                 sin(cosmo.evol.chi(i:n_z-1)),/double)
         endfor
      endfor
   end
   cosmo.const.omega lt 1.: begin
      for j=0,n_zbins-1 do begin
         g1(*,j)=2.*sinh(cosmo.evol.chi)*sinh(chi_s(j)-cosmo.evol.chi)/sinh(chi_s(j))
         for i=1,n_z-2 do begin
            g(i,j)=2.*cosmo.const.r0*$
                   int_tabulated(cosmo.evol.chi(i:n_z-1),pchi(i:n_z-1,j)*$
                                 sinh(cosmo.evol.chi(i))*$
                                 sinh(cosmo.evol.chi(i:n_z-1)-cosmo.evol.chi(i))/$
                                 sinh(cosmo.evol.chi(i:n_z-1)),/double)
         endfor
      endfor
   end
endcase
g1=g1>0.d                       ; set negative values to zero

; set g=g1 if requested
if (survey.dndztype eq 'plane') then g=g1
;if keyword_set(sheet) then g=g1
;stop
source={z_m:z_m,g:g,g1:g1,pz:pz,pchi:pchi}
return, source
end


