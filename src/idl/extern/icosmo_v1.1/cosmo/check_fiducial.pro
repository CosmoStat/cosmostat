function check_fiducial,fid,silent=silent,err_mess=err_mess

;Written by Adam Amara 15 Aug 2008
;PURPOSE: This routine checks that the enteries in the fiducial
;structure are within exceptable bounds.
;INPUTS: fiducial structure created by the routine set_fiducial
;OPTIONAL INPUTS: silent - 0: errors printed to screen, 1: no errors printed
;OUTPUTS : 1: no error detected, -1: error detected
;OPTIONAL OUTPUTS: err_mess: vector of error messages that have been collected

;*** set up initial err_message element ***
;if not keyword_set(silent) then silent=0
on_error,2
err_mess=''

;*** Check that fiducial has been input and that it is a structure ***
type=size(fid,/type)
if (type eq 0) then begin
   err_mess='ERROR: must input a fiducial structure'
    if not keyword_set(silent) then print,err_mess
    return,-1
endif 
if (type ne 8) then begin
   err_mess='ERROR: must input a fiducial must be a structure, the type input is: '+size(fid,/tname)
    if not keyword_set(silent) then print,err_mess
    return,-1
endif 


;*** check that the three structures cosmo, calc and expt are present ***
if not (tag_exist(fid,'cosmo')) then begin
   err_mess=[err_mess,'ERROR: structure cosmo is not present in fiducial structure']
   if (size(fid.cosmo,/type) ne 8) then $
      err_mess=[err_mess,'ERROR: tag cosmo in fiducial is not a structure, it is: ' $
                +size(fid.cosmo,/tname)]
endif

if not (tag_exist(fid,'calc')) then begin
   err_mess=[err_mess,'ERROR: structure calc is not present in fiducial structure']   
   if (size(fid.calc,/type) ne 8) then $
      err_mess=[err_mess,'ERROR: tag calc in fiducial is not a structure, it is: ' $
                +size(fid.calc,/tname)]
endif   

if not (tag_exist(fid,'expt')) then begin
   err_mess=[err_mess,'ERROR: structure expt is not present in fiducial structure']   
   if (size(fid.expt,/type) ne 8) then $
      err_mess=[err_mess,'ERROR: tag expt in fiducial is not a structure, it is: ' $
                +size(fid.expt,/tname)]
endif

;*** if errors have been detected then return error messages to user ***
if (n_elements(err_mess) gt 1) then begin
   err_mess=err_mess(1:*)
   if not keyword_set(silent) then print,err_mess
   return,-1
endif

;*** check cosmo structure ***
;check h:
if not (tag_exist(fid.cosmo,'h')) then begin
   err_mess=[err_mess,'ERROR: h is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.h,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: h should be in double precision instead of '+size(fid.cosmo.h,/tname)]   
   endif else begin
      if (fid.cosmo.h le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter h must be positive. It is currently set to '$
                      + strtrim(string(fid.cosmo.h),1)]
   endelse
endelse
;check omega_b:
if not (tag_exist(fid.cosmo,'omega_b')) then begin  
   err_mess=[err_mess,'ERROR: omega_b is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.omega_b,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: omega_b should be in double precision instead of '+size(fid.cosmo.omega_b,/tname)]   
   endif else begin
      if (fid.cosmo.omega_b le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_b must be positive. It is currently set to '$
                      + strtrim(string(fid.cosmo.omega_b),1)]
      if (fid.cosmo.omega_b gt 2.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_b should be roughly 0.05, it is currently set to ' $
                      + strtrim(string(fid.cosmo.omega_b),1) +'. Make sure omega_b < 2']
   endelse
endelse
;check omega_m:
if not (tag_exist(fid.cosmo,'omega_m')) then begin  
   err_mess=[err_mess,'ERROR: omega_m is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.omega_m,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: omega_m should be in double precision instead of '+size(fid.cosmo.omega_m,/tname)]   
   endif else begin
      if (fid.cosmo.omega_m le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_m must be positive. It is currently set to '$
                      + strtrim(string(fid.cosmo.omega_m),1)]
      if (fid.cosmo.omega_m lt fid.cosmo.omega_b) then $
         err_mess=[err_mess,'ERROR: Omega_m must be greater than or equal to Omega_b (Omega_m is density of all matter components, while Omega_b is the density of one component of matter - baryons)']
      if (fid.cosmo.omega_m gt 2.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_m should be roughly 0.3. It is currently set to ' $
                      + strtrim(string(fid.cosmo.omega_m),1) +'. Make sure omega_m < 2']
   endelse
endelse
;check omega_l:
if not (tag_exist(fid.cosmo,'omega_l')) then begin  
   err_mess=[err_mess,'ERROR: omega_l is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.omega_l,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: omega_l should be in double precision instead of '+size(fid.cosmo.omega_l,/tname)]   
   endif else begin
      if (fid.cosmo.omega_l le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_l must be positive. It is currently set to ' $
                   + strtrim(string(fid.cosmo.omega_l),1)]
      if (fid.cosmo.omega_l gt 2.0) then $
         err_mess=[err_mess,'ERROR: parameter omega_l should be roughly 0.7. It is currently set to ' $
                   + strtrim(string(fid.cosmo.omega_l),1) +'. Make sure omega_l < 2']     
   endelse
endelse
;check curv:
if not (tag_exist(fid.cosmo,'curv')) then begin  
   err_mess=[err_mess,'ERROR: curv is not included in the cosmo structure of fiducial structure']
endif else begin
   if ((fid.cosmo.curv eq 0) and ((fid.cosmo.omega_l + fid.cosmo.omega_m) ne 1)) then $
      err_mess=[err_mess,'ERROR: If flat Universe is set (curv=0) then the following is required: Omega_m + Omega_l = 1']
endelse
;check w0:
if not (tag_exist(fid.cosmo,'w0')) then begin  
   err_mess=[err_mess,'ERROR: w0 is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.w0,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: w0 should be in double precision instead of '+size(fid.cosmo.w0,/tname)]   
   endif else begin
      if (abs(fid.cosmo.w0) gt 3.0) then err_mess= $
         [err_mess,'ERROR: parameter w0 should be roughly -1. It is currently set to '$
          + strtrim(string(fid.cosmo.w0),1) +'. Make sure |w0| < 3']
   endelse
endelse
;check wa:
if not (tag_exist(fid.cosmo,'wa')) then begin  
   err_mess=[err_mess,'ERROR: wa is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.wa,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: wa should be in double precision instead of '+size(fid.cosmo.wa,/tname)]   
   endif else begin
      if (abs(fid.cosmo.wa) gt 2.0) then err_mess=$
         [err_mess,'ERROR: parameter wa should be roughly 0.0. It is currently set to '$
          + strtrim(string(fid.cosmo.wa),1) +'. Make sure |wa| < 2']
   endelse
endelse
;check n:
if not (tag_exist(fid.cosmo,'n')) then begin  
   err_mess=[err_mess,'ERROR: n is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.n,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: n should be in double precision instead of '+size(fid.cosmo.n,/tname)]   
   endif else begin
      if (fid.cosmo.n le 0.0) then err_mess=[err_mess,'ERROR: n should be positive. It is currently set to '$
          + strtrim(string(fid.cosmo.n),1)]
      if (fid.cosmo.n gt 3.0) then err_mess=[err_mess,'ERROR: n should be roughly 1. It is currently set to '$
          + strtrim(string(fid.cosmo.n),1)+'. Make sure 0.0 < n < 3.0']
   endelse
endelse
;check tau:
if not (tag_exist(fid.cosmo,'tau')) then begin  
   err_mess=[err_mess,'ERROR: tau is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.tau,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: tau should be in double precision instead of '+size(fid.cosmo.tau,/tname)]   
   endif else begin
      if (fid.cosmo.tau le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter tau must be positive. It is currently set to '$
          + strtrim(string(fid.cosmo.tau),1)]
      if (fid.cosmo.tau gt 1.0) then $
         err_mess=[err_mess,'ERROR: Tau should be roughly 0.09. It is currently set to '$
          + strtrim(string(fid.cosmo.tau),1)+'. Make sure Tau < 1.0']
   endelse
endelse
;check sigma8:
if not (tag_exist(fid.cosmo,'sigma8')) then begin  
   err_mess=[err_mess,'ERROR: sigma8 is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.cosmo.sigma8,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: sigma8 should be in double precision instead of '+size(fid.cosmo.sigma8,/tname)]   
   endif else begin
      if (fid.cosmo.sigma8 le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter sigma8 must be positive. It is currently set to '$
          + strtrim(string(fid.cosmo.sigma8),1)]
      if (fid.cosmo.sigma8 gt 5.0) then $
         err_mess=[err_mess,'ERROR: sigma8 should be roughly 1. It is currently set to '$
          + strtrim(string(fid.cosmo.sigma8),1)+'. Make sure sigma8 < 5.0']
   endelse
endelse
;*** Check calc structure ***
;check speed:
if not (tag_exist(fid.calc,'speed')) then begin  
   err_mess=[err_mess,'ERROR: speed is not included in the cosmo structure of fiducial structure']
endif 
;check ran_z:
if not (tag_exist(fid.calc,'ran_z')) then begin  
   err_mess=[err_mess,'ERROR: ran_z is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.calc.ran_z,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: ran_z should be in double precision instead of '+size(fid.calc.ran_z,/tname)]   
   endif else begin
      if ((where(fid.calc.ran_z lt 0.0))(0) ne -1) then $
         err_mess=[err_mess,'ERROR: ran_z must be positive. It is currently set to ['$
                      + strjoin(strtrim(string(fid.calc.ran_z)))+']']
      if (n_elements(fid.calc.ran_z) ne 2) then $
         err_mess=[err_mess,'ERROR: ran_z must a two element array (e.g. ran_z=[0.,5.]). It is currently set to '$
                      + strjoin(strtrim(string(fid.calc.ran_z)))+']']
      if (n_elements(ran_z) eq 2) then begin
         if (fid.calc.ran_z(0) ge fid.calc.ran_z(1)) then $
            err_mess=[err_mess,'ERROR: ran_z must of the form [zmin, zmax], where zmax greater than zmin. It is currently set to ['$
                         + strjoin(strtrim(string(fid.calc.ran_z)))+']']
      endif
      if (max(fid.calc.ran_z) gt 10.) then $
         err_mess=[err_mess,'ERROR: Maximum redshift (ran_z) must be smaller than z=10. It is currently set to ['$
                      + strjoin(strtrim(string(fid.calc.ran_z)))+'].']
   endelse
endelse
;check nz_fn:
if not (tag_exist(fid.calc,'nz_fn')) then begin  
   err_mess=[err_mess,'ERROR: nz_fn is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.nz_fn lt 100.0) then $
      err_mess=[err_mess,'ERROR: nz_fn must be positive and should be greater than 100. It is currently set to '$
                   + strtrim(string(fid.calc.nz_fn),1)]
   if (fid.calc.nz_fn gt 100000.0) then $
      err_mess=[err_mess,'ERROR: nz_fn is too big and will slow down computation. It is currently set to '$
                   + strtrim(string(fid.calc.nz_fn),1) +'. Make sure it is less than 100000']
   if (fid.calc.nz_fn ne round(fid.calc.nz_fn)) then $
      err_mess=[err_mess,'ERROR: nz_fn must an integer. It is currently set to '$
                   + strtrim(string(fid.calc.nz_fn),1)]
endelse
;check nz_crs:
if not (tag_exist(fid.calc,'nz_crs')) then begin  
   err_mess=[err_mess,'ERROR: nz_crs is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.nz_crs lt 20.0) then $
      err_mess=[err_mess,'ERROR: nz_crs must be positive and should be greater than 20. It is currently set to '$
                   + strtrim(string(fid.calc.nz_crs),1)]
   if (fid.calc.nz_crs gt fid.calc.nz_fn) then $
      err_mess=[err_mess,'ERROR: nz_fn must be smaller than nz_crs. They are currently set to, nz_crs= '$
                   + strtrim(string(fid.calc.nz_crs),1) +' and nz_fn= '+ strtrim(string(fid.calc.nz_fn),1) ]
   frac=fid.calc.nz_fn/fid.calc.nz_crs
   if (frac ne round(frac)) then $
      err_mess=[err_mess,'ERROR: nz_crs should be an integer of nz_fn (i.e. nz_fn/nz_crs should be a round number). They are currently set to, nz_crs= '$
                   + strtrim(string(fid.calc.nz_crs),1) +' and nz_fn= '+ strtrim(string(fid.calc.nz_fn),1) ]
endelse
;check fit_nl:
if not (tag_exist(fid.calc,'fit_nl')) then begin  
   err_mess=[err_mess,'ERROR: fit_nl is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.fit_nl ge 3.0) then $
      err_mess=[err_mess,'ERROR: fit_nl can only be set to 0,1 or 2. It is currently set to '$
                   + strtrim(string(fid.calc.fit_nl),1)]
   if (fid.calc.fit_nl ne round(fid.calc.fit_nl)) then $
      err_mess=[err_mess,'ERROR: fit_nl should be a round number. It is currently set to '$
                   + strtrim(string(fid.calc.fit_nl),1)]
endelse
;check fit_tk:
if not (tag_exist(fid.calc,'fit_tk')) then begin  
   err_mess=[err_mess,'ERROR: fit_tk is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.fit_tk ge 3.0) then $
      err_mess=[err_mess,'ERROR: fit_tk can only be set to 0,1 or 2. It is currently set to '$
                   + strtrim(string(fid.calc.fit_tk),1)]
   if (fid.calc.fit_tk ne round(fid.calc.fit_tk)) then $
      err_mess=[err_mess,'ERROR: fit_tk should be a round number. It is currently set to '$
                   + strtrim(string(fid.calc.fit_tk),1)]
endelse
;check linear:
if not (tag_exist(fid.calc,'linear')) then begin  
   err_mess=[err_mess,'ERROR: linear is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.linear gt 1.0) then $
      err_mess=[err_mess,'ERROR: linear can only be set to 0,1. It is currently set to '$
                   + strtrim(string(fid.calc.linear),1)]
   if (fid.calc.linear ne round(fid.calc.linear)) then $
      err_mess=[err_mess,'ERROR: linear should be a round number. It is currently set to '$
                   + strtrim(string(fid.calc.linear),1)]
endelse
;check k_ran:
if not (tag_exist(fid.calc,'k_ran')) then begin  
   err_mess=[err_mess,'ERROR: k_ran is not included in the cosmo structure of fiducial structure']
endif else begin
   if (size(fid.calc.k_ran,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: k_ran should be in double precision instead of '+size(fid.calc.k_ran,/tname)]   
   endif else begin
      if ((where(fid.calc.k_ran le 0.0))(0) ne -1) then $
         err_mess=[err_mess,'ERROR: k_ran must be positive. It is currently set to ['$
                      + strjoin(strtrim(string(fid.calc.k_ran)))+']']
      if (n_elements(fid.calc.k_ran) ne 2) then $
         err_mess=[err_mess,'ERROR: k_ran must a two element array (e.g. k_ran=[0.001,10.]). It is currently set to ['$
                      + strjoin(strtrim(string(fid.calc.k_ran))),']']
      if (n_elements(fid.calc.k_ran) eq 2) then begin
         if (fid.calc.k_ran(0) ge fid.calc.k_ran(1)) then $
            err_mess=[err_mess,'ERROR: k_ran must of the form [kmin, kmax], where zmax greater than zmin. It is currently set to ['$
                         + strjoin(strtrim(string(fid.calc.k_ran)))+']']
         if (fid.calc.k_ran(0) gt 0.001 ) then $
            err_mess=[err_mess,'Error: k_ran should of the form [kmin, kmax], where kmin is less than 0.001 ['$
                         + strjoin(strtrim(string(fid.calc.k_ran)))+']']
         if (fid.calc.k_ran(1) lt 100.) then $
            err_mess=[err_mess,'Error: k_ran should of the form [kmin, kmax], where kmax is greater than 100. It is currently set to ['$
                         + strjoin(strtrim(string(fid.calc.k_ran)))+']']

      endif
   endelse
endelse
;check n_k:
if not (tag_exist(fid.calc,'n_k')) then begin  
   err_mess=[err_mess,'ERROR: n_k is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.n_k lt 50.0) then $
      err_mess=[err_mess,'ERROR: n_k must be positive and should be greater than 50. It is currently set to '$
                   + strtrim(string(fid.calc.n_k),1)]
   if (fid.calc.n_k gt 1000000.) then $
      err_mess=[err_mess,'ERROR: n_k is too big and will slow down the calculation. It is currently set to n_k = '$
                   + strtrim(string(fid.calc.n_k),1)+'. Make sure n_k < 1000000']
   if (fid.calc.n_k ne round(fid.calc.n_k)) then $
      err_mess=[err_mess,'ERROR: n_k should be a round number. It is currently set to, n_k= '$
                   + strtrim(string(fid.calc.n_k),1) ]
endelse
;check l_ran:
if not (tag_exist(fid.calc,'l_ran')) then begin  
   err_mess=[err_mess,'ERROR: l_ran is not included in the cosmo structure of fiducial structure']
endif else begin
   if ((where(fid.calc.l_ran le 0.0))(0) ne -1) then $
      err_mess=[err_mess,'ERROR: l_ran must be positive. It is currently set to ['$
                   + strjoin(strtrim(string(fid.calc.l_ran)))+']']
   if (n_elements(fid.calc.l_ran) ne 2) then $
      err_mess=[err_mess,'ERROR: l_ran must a two element array (e.g. l_ran=[10.,10000.]). It is currently set to ['$
                   + strjoin(strtrim(string(fid.calc.l_ran)))+']']
   if (n_elements(l_ran) eq 2) then begin
      if (fid.calc.l_ran(0) ge fid.calc.l_ran(1)) then $
         err_mess=[err_mess,'ERROR: l_ran must of the form [lmin, lmax], where zmax greater than zmin. It is currently set to ['$
                      + strjoin(strtrim(string(fid.calc.l_ran)))+']']
   endif
endelse
;check n_l:
if not (tag_exist(fid.calc,'n_l')) then begin  
   err_mess=[err_mess,'ERROR: n_l is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.n_l lt 50.0) then $
      err_mess=[err_mess,'ERROR: n_l must be positive and should be greater than 50. It is currently set to '$
                   + strtrim(string(fid.calc.n_l),1)]
   if (fid.calc.n_l gt 1000000.) then $
      err_mess=[err_mess,'ERROR: n_l is too big and will slow down the calculation. It is currently set to n_l = '$
                   + strtrim(string(fid.calc.n_l),1)+'. Make sure n_l < 1000000']
   if (fid.calc.n_l ne round(fid.calc.n_l)) then $
      err_mess=[err_mess,'ERROR: n_l should be a round number. It is currently set to, n_l= '$
                   + strtrim(string(fid.calc.n_l),1) ]
endelse
;check n_lbin:
if not (tag_exist(fid.calc,'n_lbin')) then begin  
   err_mess=[err_mess,'ERROR: n_lbin is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.n_lbin le 0.0) then $
      err_mess=[err_mess,'ERROR: n_lbin must be positive. It is currently set to '$
                   + strtrim(string(fid.calc.n_lbin),1)]
   if (fid.calc.n_lbin gt fid.calc.n_l) then $
      err_mess=[err_mess,'ERROR: n_lbin should be smaller than b_l. They are currently set to n_lbin = '$
                   + strtrim(string(fid.calc.n_lbin),1)+' and n_l =  ' + strtrim(string(fid.calc.n_l),1)]
   if (fid.calc.n_lbin ne round(fid.calc.n_lbin)) then $
      err_mess=[err_mess,'ERROR: n_lbin should be a round number. It is currently set to, n_lbin= '$
                   + strtrim(string(fid.calc.n_lbin),1) ]
endelse
;check delta:
if not (tag_exist(fid.calc,'delta')) then begin  
   err_mess=[err_mess,'ERROR: delta is not included in the calc structure of fiducial structure']
endif else begin
   if (size(fid.calc.delta,/type) ne 5) then begin
      err_mess=[err_mess,'ERROR: delta should be in double precision instead of '+size(fid.calc.delta,/tname)]   
   endif else begin
      if (fid.calc.delta le 0.0) then $
         err_mess=[err_mess,'ERROR: parameter delta must be positive. It is currently set to '$
                      + strtrim(string(fid.calc.delta),1)]
      if (fid.calc.delta gt 1.0) then $
         err_mess=[err_mess,'ERROR: delta should be roughly 0.003. It is currently set to '$
                      + strtrim(string(fid.calc.delta),1)+'. Make sure delta < 1.0']
   endelse
endelse
;check verbose:
if not (tag_exist(fid.calc,'verbose')) then begin  
   err_mess=[err_mess,'ERROR: verbose is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.verbose gt 2.0) then $
      err_mess=[err_mess,'ERROR: verbose can only be set to 0, 1 or 2. It is currently set to '$
                   + strtrim(string(fid.calc.verbose),1)]
   if (fid.calc.verbose ne round(fid.calc.verbose)) then $
      err_mess=[err_mess,'ERROR: verbose should be a round number. It is currently set to '$
                   + strtrim(string(fid.calc.verbose),1)]
endelse
;check err:
if not (tag_exist(fid.calc,'err')) then begin  
   err_mess=[err_mess,'ERROR: err is not included in the cosmo structure of fiducial structure']
endif else begin
   if (fid.calc.err gt 2.0) then $
      err_mess=[err_mess,'ERROR: err can only be set to 0, 1 or 2. It is currently set to '$
                   + strtrim(string(fid.calc.err),1)]
   if (fid.calc.err ne round(fid.calc.err)) then $
      err_mess=[err_mess,'ERROR: err should be a round number. It is currently set to '$
                   + strtrim(string(fid.calc.err),1)]
endelse

;*** Check expt structure ***
svs=tag_names(fid.expt)
for i=0,n_elements(svs)-1 do begin
;check n_zbin
   if not (tag_exist(fid.expt.(i),'n_zbin')) then begin  
      err_mess=[err_mess,'ERROR: n_zbin is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (fid.expt.(i).n_zbin le 0.0) then $
         err_mess=[err_mess,'ERROR: n_zbin must be positive. It is currently set to '$
                   + strtrim(string(fid.expt.(i).n_zbin),1)+' in '+svs(i)]
      if (fid.expt.(i).n_zbin gt 200) then $
         err_mess=[err_mess,'ERROR: n_zbin is too big. It is currently set to '$
                   + strtrim(string(fid.expt.(i).n_zbin),1)+' in expt.'+svs(i)+'. Make sure it is less than 200']
      if (round(fid.expt.(i).n_zbin) ne  fid.expt.(i).n_zbin) then $
         err_mess=[err_mess,'ERROR: n_zbin should be a round number. It is currently set to '$
                   + strtrim(string(fid.expt.(i).n_zbin),1)+' in expt.'+svs(i)]
   endelse
;check zerror:
   if not (tag_exist(fid.expt.(i),'zerror')) then begin  
      err_mess=[err_mess,'ERROR: zerror is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).zerror,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: zerror should be in double precision instead of '+size(fid.expt.(i).zerror,/tname)]   
      endif else begin
         if (fid.expt.(i).zerror lt 0.0) then $
            err_mess=[err_mess,'ERROR: zerror must be positive. It is currently set to '$
                      + strtrim(string(fid.expt.(i).zerror),1)+' in  expt.'+svs(i)]
         if (fid.expt.(i).zerror gt 1.0) then $
            err_mess=[err_mess,'ERROR: zerror should be ~ 0.05. It is currently set to '$
                      + strtrim(string(fid.expt.(i).zerror),1)+' in expt.'+svs(i)+'. Make sure zerror is less than 1']
      endelse
   endelse
;check z_med:
   if not (tag_exist(fid.expt.(i),'z_med')) then begin  
      err_mess=[err_mess,'ERROR: z_med is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).z_med,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: z_med should be in double precision instead of '+size(fid.expt.(i).z_med,/tname)]   
      endif else begin
         if (fid.expt.(i).z_med lt 0.0) then $
            err_mess=[err_mess,'ERROR: z_med must be positive (or zero if not set). It is currently set to '$
                      + strtrim(string(fid.expt.(i).z_med),1)+' in expt.'+svs(i)]
         if (fid.expt.(i).z_med gt max(fid.calc.ran_z)/2.) then $
            err_mess=[err_mess,'Warning: z_med might be too large compared with ran_z. These are currently set to, z_med ='$
                      + strtrim(string(fid.expt.(i).z_med),1) +', ran_z = ' $
                      +strjoin(strtrim(string(fid.calc.ran_z)))+$
                      '. Make sure maximum range in ran_z is at least 2 * z_med']
      endelse
   endelse
;check ng:
   if not (tag_exist(fid.expt.(i),'ng')) then begin  
      err_mess=[err_mess,'ERROR: ng is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).ng,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: ng should be in double precision instead of '+size(fid.expt.(i).ng,/tname)]   
      endif else begin
         num_ng=n_elements(fid.expt.(i).ng)
         if not ((num_ng eq 1) or (num_ng eq fid.expt.(i).n_zbin)) then $
            err_mess=[err_mess,'ERROR: number of elements in n_g is wrong. It should either be 1 (i.e. total number of galaxies) or n_zbin (i.e. number of galaxies per z bin']

         if ((where(fid.expt.(i).ng lt 0.0))(0) ne -1) then $
            err_mess=[err_mess,'ERROR: ng must be positive. It is currently set to '$
                      + strjoin(strtrim(string(fid.expt.(i).ng)))+' in expt.'+svs(i)]
         if (total(fid.expt.(i).ng) le 0.0) then $
            err_mess=[err_mess,'ERROR: sum of ng must be positive. It is currently '$
                      + strtrim(string(total(fid.expt.(i).ng)),1)+' in expt.'+svs(i)]
         if (max(fid.expt.(i).ng) gt 10000.0) then $
            err_mess=[err_mess,'Error: elements of ng should be roughly 1 -> 100. it currently set to, ng = [' $
                      + strtrim(string(fid.expt.(i).ng),1) +']. Make sure all elements are smaller than 10000.0']
      endelse
   endelse
;check a_survey:
   if not (tag_exist(fid.expt.(i),'a_survey')) then begin  
      err_mess=[err_mess,'ERROR: a_survey is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).a_survey,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: a_survey should be in double precision instead of '+size(fid.expt.(i).a_survey,/tname)]   
      endif else begin
         if (fid.expt.(i).a_survey le 0.0) then $
            err_mess=[err_mess,'ERROR: a_survey must be positive. It is currently set to '$
                      + strtrim(string(fid.expt.(i).a_survey),1)+' in '+svs(i)]
         if (fid.expt.(i).a_survey gt 41253.) then $
            err_mess=[err_mess,'ERROR: a_survey must be less than ~41253 (square deg in full sky). It is currently set to '$
                      + strtrim(string(fid.expt.(i).a_survey),1)+' in '+svs(i)]
      endelse
   endelse
 ;check eff:
   if not (tag_exist(fid.expt.(i),'eff')) then begin  
      err_mess=[err_mess,'ERROR: eff is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).eff,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: eff should be in double precision instead of '+size(fid.expt.(i).eff,/tname)]   
      endif else begin
         if (fid.expt.(i).eff le 0.0) then $
            err_mess=[err_mess,'ERROR: eff must be positive. It is currently set to '$
                      + strtrim(string(fid.expt.(i).eff),1)+' in '+svs(i)]
         if (fid.expt.(i).eff gt 1.) then $
            err_mess=[err_mess,'ERROR: eff must be less than 1. It is currently set to '$
                      + strtrim(string(fid.expt.(i).eff),1)+' in '+svs(i)]
      endelse
   endelse
 ;check sig_int:
   if not (tag_exist(fid.expt.(i),'sig_int')) then begin  
      err_mess=[err_mess,'ERROR: sig_int is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).sig_int,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: sig_int should be in double precision instead of '+size(fid.expt.(i).sig_int,/tname)]   
      endif else begin
         if (fid.expt.(i).sig_int le 0.0) then begin
            if (fid.expt.(i).sig_int ne -1) then $
               err_mess=[err_mess,'ERROR: sig_int must be positive (or -1 if not set). It is currently set to '$
                         + strtrim(string(fid.expt.(i).sig_int),1)+' in '+svs(i)]
         endif
         if (fid.expt.(i).sig_int gt 10.) then $
            err_mess=[err_mess,'ERROR: sig_int should be roughly 0.25. It is currently set to '$
                      + strtrim(string(fid.expt.(i).sig_int),1)+' in '+svs(i)+'. Make sure sig_int is less than 10']
      endelse
   endelse
;check dndztype:
   if not (tag_exist(fid.expt.(i),'dndztype')) then begin  
      err_mess=[err_mess,'ERROR: dndztype is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).dndztype,/type) ne 7) then begin
         err_mess=[err_mess,'ERROR: dndztype should be a string instead of '+size(fid.expt.(i).dndztype,/tname)]   
      endif else begin
         poss_types=['smail','plane','hist']
         if ((where(strlowcase(poss_types) eq fid.expt.(i).dndztype))(0) eq -1) then $
            err_mess=[err_mess,'ERROR: dndztype must either: smail or plane or hist it is currently set to'$
                      + strtrim(string(fid.expt.(i).dndztype),1)+' in '+svs(i)]
      endelse
   endelse

;check biastype:
   if not (tag_exist(fid.expt.(i),'biastype')) then begin  
      err_mess=[err_mess,'ERROR: biastype is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).biastype,/type) ne 7) then begin
         err_mess=[err_mess,'ERROR: biastype should be a string instead of '+size(fid.expt.(i).biastype,/tname)]   
      endif else begin
         poss_types=['bias1']
         if ((where(strlowcase(poss_types) eq fid.expt.(i).biastype))(0) eq -1) then $
            err_mess=[err_mess,'ERROR: biastype must either: smail or plane or hist it is currently set to'$
                      + strtrim(string(fid.expt.(i).biastype),1)+' in '+svs(i)]
      endelse
   endelse
 ;check ns:
   if not (tag_exist(fid.expt.(i),'ns')) then begin  
      err_mess=[err_mess,'ERROR: ns is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).ns,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: ns should be in double precision instead of '+size(fid.expt.(i).ns,/tname)]   
      endif else begin
         if (fid.expt.(i).ns le 0.0) then begin
            if (fid.expt.(i).ns ne -1) then $
               err_mess=[err_mess,'ERROR: ns must be positive (or -1 if not set). It is currently set to '$
                         + strtrim(string(fid.expt.(i).ns),1)+' in '+svs(i)]
         endif
         if (fid.expt.(i).ns gt 10.) then $
            err_mess=[err_mess,'ERROR: ns should be roughly 0.03. It is currently set to '$
                      + strtrim(string(fid.expt.(i).ns),1)+' in '+svs(i)+'. Make sure ns is less than 10']
      endelse
   endelse
 ;check sigmam:
   if not (tag_exist(fid.expt.(i),'sigmam')) then begin  
      err_mess=[err_mess,'ERROR: sigmam is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).sigmam,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: sigmam should be in double precision instead of '+size(fid.expt.(i).sigmam,/tname)]   
      endif else begin
         if (fid.expt.(i).sigmam le 0.0) then begin
            if (fid.expt.(i).sigmam ne -1) then $
               err_mess=[err_mess,'ERROR: sigmam must be positive (or -1 if not set). It is currently set to '$
                         + strtrim(string(fid.expt.(i).sigmam),1)+' in '+svs(i)]
         endif
         if (fid.expt.(i).sigmam gt 100.) then $
            err_mess=[err_mess,'ERROR: sigmam should be roughly 0.2. It is currently set to '$
                      + strtrim(string(fid.expt.(i).sigmam),1)+' in '+svs(i)+'. Make sure sigmam is less than 100']
      endelse
   endelse
 ;check delm:
   if not (tag_exist(fid.expt.(i),'delm')) then begin  
      err_mess=[err_mess,'ERROR: delm is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).delm,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: delm should be in double precision instead of '+size(fid.expt.(i).delm,/tname)]   
      endif else begin
         if (fid.expt.(i).delm le 0.0) then begin
            if (fid.expt.(i).delm ne -1) then $
               err_mess=[err_mess,'ERROR: delm must be positive (or -1 if not set). It is currently set to '$
                         + strtrim(string(fid.expt.(i).delm),1)+' in '+svs(i)]
         endif
         if (fid.expt.(i).delm gt 10.) then $
            err_mess=[err_mess,'ERROR: delm should be roughly 0.02. It is currently set to '$
                      + strtrim(string(fid.expt.(i).delm),1)+' in '+svs(i)+'. Make sure delm is less than 10']
      endelse
   endelse
;check sne_zran:
   if not (tag_exist(fid.expt.(i),'sne_zran')) then begin  
      err_mess=[err_mess,'ERROR: sne_zran is not included in the cosmo structure of fiducial structure'+' in '+svs(i)]
   endif else begin
      if (size(fid.expt.(i).sne_zran,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: sne_zran should be in double precision instead of '+size(fid.expt.(i).sne_zran,/tname)]   
      endif else begin
         if (fid.expt.(i).sne_zran(0) ne -1) then begin
            if ((where(fid.expt.(i).sne_zran le 0.0))(0) ne -1) then $
               err_mess=[err_mess,'ERROR: sne_zran must be positive. It is currently set to ['$
                         + strjoin(strtrim(string(fid.expt.(i).sne_zran)))+']']
            if (n_elements(fid.expt.(i).sne_zran) ne 2) then $
               err_mess=[err_mess,'ERROR: sne_zran must a two element array (e.g. sne_zran=[0.,5.]). It is currently set to '$
                         + strjoin(strtrim(string(fid.expt.(i).sne_zran)))+']']
            if (n_elements(fid.expt.(i).sne_zran) eq 2) then begin
               if (fid.expt.(i).sne_zran(0) gt fid.expt.(i).sne_zran(1)) then $
                  err_mess=[err_mess,'ERROR: sne_zran must of the form [zmin, zmax], where zmax greater than zmin. It is currently set to ['$
                            + strjoin(strtrim(string(fid.expt.(i).sne_zran)))+']']
            endif
               if (max(fid.expt.(i).sne_zran) gt max(fid.calc.ran_z)) then $
                  err_mess=[err_mess,'ERROR: Maximum redshift (sne_zran) must be smaller than the calc.ran_z. It is currently set to ['$
                            + strjoin(strtrim(string(fid.expt.(i).sne_zran)))+']']
         endif
      endelse
   endelse
;check name:
   if not (tag_exist(fid.expt.(i),'name')) then begin  
      err_mess=[err_mess,'ERROR: name is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).name,/type) ne 7) then begin
         err_mess=[err_mess,'ERROR: name should be a string instead of '+size(fid.expt.(i).name,/tname)]   
      endif 
   endelse

;check probes:
   if not (tag_exist(fid.expt.(i),'probes')) then begin  
      err_mess=[err_mess,'ERROR: probes is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).probes,/type) ne 7) then begin
         err_mess=[err_mess,'ERROR: probes should be a string instead of '+size(fid.expt.(i).probes,/tname)]   
      endif else begin
         poss_types=['lens','bao','sne']
         for j=0,n_elements(fid.expt.(i).probes)-1 do begin
            if ((where(strlowcase(poss_types) eq fid.expt.(i).probes(j)))(0) eq -1) then $
               err_mess=[err_mess,'ERROR: probes must either: smail or plane or hist it is currently set to'$
                         + strtrim(string(fid.expt.(i).probes(j)),1)+' in '+svs(i)]
         endfor
      endelse
   endelse
;check dndzp:
   if not (tag_exist(fid.expt.(i),'dndzp')) then begin  
      err_mess=[err_mess,'ERROR: dndzp is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).dndzp,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: dndzp should be in double precision instead of '+size(fid.expt.(i).dndzp,/tname)]   
      endif else begin
         num_dndzp=n_elements(fid.expt.(i).dndzp)
         if (strlowcase(fid.expt.(i).dndztype) eq 'smail') then begin
            if not (num_dndzp eq 2) then $
               err_mess=[err_mess,'ERROR: number of elements in dndzp is wrong. For smail et al this should be 2']   
         endif         
         if ((where(fid.expt.(i).dndzp lt 0.0))(0) ne -1) then begin
            if (fid.expt.(i).dndzp(0) ne -1) then $
            err_mess=[err_mess,'ERROR: dndzp must be positive. It is currently set to '$
                      + strjoin(strtrim(string(fid.expt.(i).dndzp)))+' in expt.'+svs(i)]
         endif
      endelse
   endelse
;check dndzz:
   if not (tag_exist(fid.expt.(i),'dndzz')) then begin  
      err_mess=[err_mess,'ERROR: dndzz is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
   endif else begin
      if (size(fid.expt.(i).dndzz,/type) ne 5) then begin
         err_mess=[err_mess,'ERROR: dndzz should be in double precision instead of '+size(fid.expt.(i).dndzz,/tname)]   
      endif else begin
         num_dndzz=n_elements(fid.expt.(i).dndzz)
         if ((where(fid.expt.(i).dndzz lt 0.0))(0) ne -1) then begin
            if (fid.expt.(i).dndzz(0) ne -1) then $
               err_mess=[err_mess,'ERROR: dndzz must be positive. It is currently set to '$
                         + strjoin(strtrim(string(fid.expt.(i).dndzz)))+' in expt.'+svs(i)]
         endif
      endelse
   endelse

endfor

;*** if errors have been detected then return error message to user ***
if (n_elements(err_mess) gt 1) then begin
   err_mess=err_mess(1:*)
   if not keyword_set(silent) then print,transpose(err_mess)
   return,-1
endif else begin
   err_mess=0
   return,1
endelse



;DNDZP
;DNDZZ



;;check ***:
;if not (tag_exist(fid.expt.(i),'***')) then begin  
;   err_mess=[err_mess,'ERROR: *** is not included in the '+svs(i)+'structure in expt structure of fiducial structure'];
;endif else begin
;   if (type ne 5) then begin
;      err_mess=[err_mess,'ERROR: *** should be in double precision instead of '+size(fid.expt.(i).***,/tname)]   
;   endif else begin
;      if (fid.expt.(i).*** le 0.0) then $
;         err_mess=[err_mess,'ERROR: *** must be positive. It is currently set to '$
;                   + strtrim(string(fid.expt.(i).***),1)+' in '+svs(i)]
;   endelse
;endelse

;check ***:
;if not (tag_exist(fid.expt.(i),'***')) then begin  
;   err_mess=[err_mess,'ERROR: *** is not included in the '+svs(i)+'structure in expt structure of fiducial structure']
;endif else begin
;   if (fid.calc.*** le 0.0) then $
;      err_mess=[err_mess,'ERROR: *** must be positive. It is currently set to '$
;                + strtrim(string(fid.calc.***),1)]
;endelse


end


