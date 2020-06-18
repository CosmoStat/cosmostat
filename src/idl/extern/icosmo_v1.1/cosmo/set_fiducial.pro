function set_fiducial,experiment,cosmo_in=cosmo_in,expt_in=expt_in,calc_in=calc_in,silent=silent

; Sep 08 modified by AnR, fixed bug in choice of experiments
; Aug 08 modified by AA, header updated
; Aug 08 modified by AA, hardwire structure elements as double precision
; Aug 08 modified by TK to include SNe options
; Aug 08 modified by AnR, fixed bug: comso_in instead of cosmo_in in
; some places
; Aug 08 modified by AA suppress alerts and warnings (using !quiet=1)
; Aug 08 modified by AA to tidy up calc options
; Jul 08 modified by AA to include user input structures
; Jul 08 modified by AA to extend the survey structure for each probe
; Jul 08 modified by AA to include set_keyword, and rearrange structure
; Written by Adam Amara 15th May 2008
; ***************// HELP BEGIN //**************
; PURPOSE: This routine is designed to set the fiducial parameter sets
;          that are used throughout the iCosmo package. The user is
;          able to set all the options here. (Note: At the moment only 
;          sigma8 normalisation is used.)
;
; OPTIONAL INPUT: 
;          experiment: Name of one of the preset
;                      experiments. Options are:
;                 *** NEED TO ADD ***
;          cosmo_in: input cosmological parameter, the options are: 
;              H: hubble constant 
;              Omega_B: baryon density 
;              Omega_m: matter density 
;              Omega_L: dark energy density 
;              w0, wa : equation of state parameters - w(a)=w0+(1-a)wa),  
;              n: spectral index 
;              Tau: optical depth 
;              Sigma8: power spectrum normalisation  
;              Curv: curvature: 0- flat & 1-curved 
;
;          calc_in: Input calculation parameters. Options are:
;              fit_nl: non-linear correction 
;                      (0: P&D, 1: Ma et al, 2:Smith et al).), 
;              fit_TK: Transfer function 
;                      (0: E&H w/o wiggles, 1: E&H w wiggles, 2: BBKS) 
;              verbose:prited comments: 
;                      0:nothing is printed, 1: some useful comments
;                      printed, 2: many comments are printed, 
;              delta:  Step size in Fisher matrix calculations.
;              n_l:    Number of points in ell range.
;              speed:  Choice of either speed or accuracy. 
;                      (0:accuarate, 1: fast)
;              linear: Use the linear power spectrum when set
;              k_ran:  Range of k values for the integration [h Mpc^-1]
;              n_k:    Number of k values for integration (default=200)
;              n_lbin: Number of l-bins (default: 20)
;              nz_fn:  number of points in creating an array of FINE z values
;              nz_crs: number of points in creating an array of COARSE z
;                      values. This should be a smaller than (or equal
;                      to) nz_fn and nz_fn/nz_crs should be a round
;                      number 
;              ran_z:  Range of z values (default: [0,5])
;              l_ran:  lensing ell range
;              err:    Error handling keyword. 0: stops where error
;                      occured, stops and returns to main level,2: stops
;                      and returns to the routine that called the
;                      routine that caused the error 
;        expt_in: Input experiment properties of a survey. Options
;                 (note that sv1 can be changed to sv2 or sv3):
;              sv1_N_ZBIN:  The number of redshift slices (lensing tomography)
;              sv1_ZERROR:  z error, sigma(z) = gerror*(1+z)
;              sv1_Z_MED:   Median redshift of all galaxies
;              sv1_NG:      number density of galaxies [per arc min^2] 
;                           (scalr or a vector (where the vector
;                           version has the same number of elements at
;                           there are tomography bin) 
;              sv1_A_SURVEY:area of survey [sq. deg.]
;              sv1_EFF:     masking efficiency
;              sv1_SIG_INT: sigma_e1/G from RRG I from intrinsic+noise
;                           <gamma_1^2>^.5=<gamma_2^2>^.5  
;              sv1_DNDZTYPE:form of galaxy distribution for lensing 
;              sv1_DNDZP:   parameters of galaxy distribution 
;                           (e.g. for Smail et al. p(0) = alpha and (1) = beta)
;              sv1_DNDZZ:   redshift information for gal dist
;              sv1_BIASTYPE:bias type 
;              sv1_NS:      number of Sne [# Sne/arcmin^2]
;              sv1_SIGMAM:  intrisic scatter in Sne mag [magnitude]
;              sv1_DELM:    magnitude precision [magnitude]
;              sv1_SNE_ZRAN:Sne redshift range
;              sv1_NAME:    name of the survey
;              sv1_PROBES:  possible probes measured by the survey
;
; OUTPUT: strucutre containing all the relevant parameters used in the
;         iCosmo package
;
; ----
; Example 1: 
; > fid=set_fiducial('generic',cosmo_in={omega_m:.25d},calc_in={fit_tk:2},expt_in{sv1_n_zbin:2,sv1_zerror:0.02d})
;
; This will set the fiducial parameters to:
; omega_m=0.25, tranfer function set to BBKS, 2 redshift bins are used
; and the photometric redshift errors are 0.02. All other parameters are
; set to the default values.

; Example 2: 
; > cosm = {h:0.8,n:1.1,curv:0}
; > calc = {nz_fn:1000,fit_nl:2,ran_z:[0.0d,7.0d]}
; > expt = {sv1_a_sruvey:10000.0d,sv1_z_med=0.8d,sv2_ng:8.0d,sv3_sig_int:0.2d}
; > fid=set_fiducial('generic',cosmo_in=cosm,calc_in=calc,expt_in=expt)
;
; In this example the three input structures are setup before
; set_fiducial.pro is called. 
; Here the cosmology parameters are set to:
;      h = 0.8, n = 1.1 and a flat model is used
; the calculation parameters are set to:
;      nz_fn = 1000, Peacock and Dodds is used for the non-linear
;      correction and the redshift range is set to [0,7]
; the experiment 'generic' is used with the following survey
; modification:
;      area and median redshift of sv1 is set to 10000 deg^2 and 0.8
;      number of galaxies in sv2 is set to 8 gal/arcmin^2
;      intrinsic shear in sv3 set to 0.2 
; ----
; ***************// HELP END //**************


;***************** Begin: set up Cosmology Structure ********************
if (keyword_set(cosmo_in) or keyword_set(expt_in) or keyword_set(calc_in)) then begin
   check=check_fidinput(cosmo_in=cosmo_in,expt_in=expt_in,calc_in=calc_in)
   if (check ne 0) then return, -1   
endif 

;***************** end: set up Cosmology Structure *********************

;***************** Begin: set up Cosmology Structure *********************  
; Hubble constant: Ho=h*100 km/s/Mpc:  
  if not tag_check(cosmo_in,'h',val=h) then h = 0.7d 
; Baryon density:
  if not tag_check(cosmo_in,'omega_b',val=omega_b) then  omega_b = 0.045d               
; Total matter density (DM + B):
  if not tag_check(cosmo_in,'omega_m',val=omega_m) then  omega_m = 0.3d                
; Dark energy density:
  if not tag_check(cosmo_in,'omega_l',val=omega_l) then  omega_l = 0.7d                
; w parametrised as: w(a)=w0+wa(1-a):
  if not tag_check(cosmo_in,'w0',val=w0) then  w0 = -0.95d                    
  if not tag_check(cosmo_in,'wa',val=wa) then  wa = 0.0d                     
; Spectral index:
  if not tag_check(cosmo_in,'n',val=n) then  n = 1.0d                       
; Optical depth (only for CMB norm.):
  if not tag_check(cosmo_in,'tau',val=tau) then  tau = 0.09d                   
; sigma(8Mpc) for normalisation:
  if not tag_check(cosmo_in,'sigma8',val=sigma8) then  sigma8 = 0.8d                  
; Set curved or flat models:
  if not tag_check(cosmo_in,'curv',val=curv) then curv = 1b 

;Check consistency of total density and curvature:
  if ((curv eq 0) and (omega_m+omega_l ne 1.)) then begin
;     stop
     if (tag_check(cosmo_in,'omega_m') and tag_check(cosmo_in,'omega_l')) then begin
        print, 'You cannot set the model to be flat and set omega_m+omega_l not equal to 1'
        return,'Flatness Error'
     endif 
     if tag_check(cosmo_in,'omega_m') then omega_l=1.d -omega_m
     if tag_check(cosmo_in,'omega_l') then omega_m=1.d -omega_l
  endif

;put fiducial values into a structure:  
cosmo={h:double(h),omega_b:double(omega_b),omega_m:double(omega_m),omega_l:double(omega_l),w0:double(w0),wa:double(wa),n:double(n),tau:double(tau),sigma8:double(sigma8),curv:long(curv)}  
;***************** End: Set Up Cosmology Structure *********************  


;***************** Begin: Set Up Calc Structure *********************  
; *** Speed option ***
; speed: Choice of either speed or accuracy. 
; Used by: ???
if not tag_check(calc_in,'speed',val=speed) then speed=1b
; *** redshift options: ***
; ran_z: Range of z values (default: [0,5]). 
; Used by: mk_evol, mk_fisher_lens
if not tag_check(calc_in,'ran_z',val=ran_z) then ran_z=[0.d,5.d]
;nz_fn: number of points in creating an array of FINE z values
;used by: mk_cosmo and mk_survey
  if not tag_check(calc_in,'nz_fn',val=nz_fn) then $
     if (speed) then nz_fn=400L else nz_fn=800L
;nz_crs: number of points in creating an array of COARSE z
;values. This should be a smaller than (or equal to) nz_fn and
;nz_fn/nz_crs should be a round number. 
;(e.g. nz_fn=100, nz_crs=20 - ok; nz_fn=100, nz_crs=21 - not ok)
;used by: mk_pk
  if not tag_check(calc_in,'nz_crs',val=nz_crs) then $
     if (speed) then nz_crs=100L else nz_crs=200L

; *** 3D powerspectrum options: ***
; fit_nl: non-linear correction (0: P&D, 1: Ma et al, 2:Smith et al). 
; Used by: mk_cl_tomo, mk_del2
  if not tag_check(calc_in,'fit_nl',val=fit_nl) then fit_nl=2b
; fit_tk: Transfer function (0: E&H w/o wiggles, 1: E&H w wiggles, 2: BBKS). 
; Used by: mk_cl_tomo, del2_lin
  if not tag_check(calc_in,'fit_tk',val=fit_tk) then fit_tk=0b
; linear: compute linear power spectrum when set. 
; Used by: mk_cl_tomo
if not tag_check(calc_in,'linear',val=linear) then linear=0
; k_ran: Range of k values for the integration [h Mpc^-1]
; Used by: & ???
if not tag_check(calc_in,'k_ran',val=k_ran) then k_ran=[.001d,1000.d]  
; n_k: Number of k values for integration (default=200). 
; Used by: mk_cl_tomo, mk_del2
if not tag_check(calc_in,'n_k',val=n_k) then $
   if (speed) then n_k=200L else n_k=400L

; *** lensing powespectrum options: ***
; l_ran: ell range.
; Used by: mk_cl_tomo, mk_cl_cov, mk_fisher_lens
  if not tag_check(calc_in,'l_ran',val=l_ran) then l_ran=[10.0d,2.0d4]
; n_l: Number of points in ell range. 
; Used by: mk_cl_tomo, mk_fisher_lens
  if not tag_check(calc_in,'n_l',val=n_l) then $
     if (speed) then n_l=200L else n_l=400L
; n_lbin: Number of l-bins (default: 20)
; Used by: mk_cl_cov_tomo
if not tag_check(calc_in,'n_lbin',val=n_lbin) then n_lbin=40L

; *** fisher matrix calculation options: ***
; delta: Step size in Fisher matrix calculations. 
; Used by: mk_cldp_4p 
  if not tag_check(calc_in,'delta',val=delta) then delta=0.003d

; *** other options: ***
; verbose: Verbose setting: 0: nothing is printed, 1: some output is
; printed (i.e usual user mode) and 2: lots of information is printed
; (i.e. developer mode). 
  if not tag_check(calc_in,'verbose',val=verbose) then verbose=1b
;err: Error handling keyword. 0: stops where error occured, stops and
;returns to main level,2: stops and returns to the routine that called
;the routine that caused the error.
  if not tag_check(calc_in,'err',val=err) then err=1b

calc= {fit_nl:long(fit_nl),fit_tk:long(fit_tk),verbose:long(verbose),delta:double(delta),n_l:n_l,l_ran:double(l_ran),speed:long(speed),linear:long(linear),k_ran:double(k_ran),n_k:long(n_k),n_lbin:long(n_lbin),nz_fn:long(nz_fn),nz_crs:long(nz_crs),ran_z:double(ran_z),err:long(err)}
;***************** End: Set Up Calc Structure *********************  

;***************** Begin: Set Up Survey Structure *********************  
;*** Set base experiment ***
;*** At the moment any experiment can have upto three surveys ***


if not keyword_set(experiment) then experiment = 'generic'

case strlowcase(experiment) of 
   'generic': expt=expt_generic(expt_in) ; generic survey setup
   'des': expt=expt_des(expt_in)         
   'essence': expt=expt_essence(expt_in)
   'lsst': expt=expt_lsst(expt_in)
   'panstarrs1': expt=expt_panstarrs1(expt_in)
   'snap':expt=expt_snap(expt_in)
   'snls': expt=expt_snls(expt_in)
   'boss': expt=expt_boss(expt_in)
   'wfmos': expt=expt_wfmos(expt_in)
   'wigglez':expt=expt_wigglez(expt_in)
   'cosmos':expt=expt_cosmos(expt_in)
   'slitless': expt=expt_slitless(expt_in)
   else: begin
      print,'The survey you requested '+experiment+'is not supported.'
      print,'Please select one of the following surveys:'
      print,'generic, DES, ESSENCE, LSST, PanSTARRS1, SNAP, SNLS,BOSS, WFMOS,Wigglez'
      return,{error:-1}
   end
endcase
;***************** End: Set Up Survey Structure *********************  

; -- start: Construct the keywords structure --
fiducial={cosmo:cosmo,expt:expt,calc:calc}
; -- end: Construct the keywords structure --

; -- start: check structure --
temp=check_fiducial(fiducial,err_mess=err_mess,silent=silent)
if (temp ne 1) then begin
   return,{error:-1}
endif
; -- end: check structure --


return,fiducial

end




