;function mk_dcldp_4p,fiducial,sv,step_param,keywords=keywords
function mk_dcldp_4p,fid,sv,step_param

; Jul 08 - Modified by AA compatible with new set_fiducial & mk_cosmo
; Written by Adam Amara 10 May 2006.
; PURPOSE: This routines calculates the 1st derivative of a function,
; with respect to one of the cosmology parameters using 4 points.
; INPUTS: fiducial: Structure with the cosmological model (created
;                   using set_fiducial)
;         sv: survey parameters structure (created using mk_survey)
;         step_param: name of parameter which will be used to
;                   step away from the fiducial values. 
; OPTIONAL INPUTS: keywords (as set by the routine keywords.pro). 
;                  The following keywords are used:'
;                      delta: step size (default: delta=0.003). For
;                         parameters with values greater than delta, the
;                         step in that parameter  (p) will be
;                         delta_p=delta*p, otherwise delta_p= delta  
;                      printit: verbose mode (default: printit=0b)
;                      ** keywords used by mk_evol,mk_csource and
;                      mk_source will be passed to these routines **
; OUTPUT: returns derivatives. 
; ----
; Example 1: dcldp_temp=mk_dcldp_4p(fiducial,sv,'h')
; Calculates dcl/dh (derivative with respect to the hubble
; parameter) using the fiducial cosmology (structure created by
; set_fiducial) and survey sv (created by mk_survey). Optional
; keywords are set to inbuilt default values. 
;
; Example 2: dcldp_temp=mk_dcldp_4p(fiducial,sv,'sigma8',keywords=keywords)
; Calculates dcl/dsigma8 (derivative with respect to the sigma8, the
; power spectrum normalisation) using the fiducial cosmology
; (structure created by  set_fiducial) and survey sv (created by
; mk_survey). Optional keywords are set by the keywords structure
; created by the set_keywords routine. Typically sigma8 of order
; unity, in which case the step size will be: delta*sigma8 
; 
; Example 3: dcldp_temp=mk_dcldp_4p(fiducial,sv,'wa',keywords={delta:0.01})
; Calculates dcl/wa (derivative with respect to the wa) using the
; fiducial cosmology (structure created by  set_fiducial) and survey
; sv (created by mk_survey). The optional keyword delta is set of
; delta=0.01. If wa=0 then the step_size will be equal to delta.
; ----


if not tag_exist(fid.cosmo,step_param,index=ind) then begin
   print,'The parameter you are trying to differentiate with respect to is not one the main cosmology parameters.'
   stop
   return,'Error: Parameter not used'
endif
if (fid.calc.verbose ge 1) then print,'varying: '+step_param 
;if not tag_check(keywords,'delta',val=delta) then delta=0.003d
delta=fid.calc.delta

if (abs(fid.cosmo.(ind)) le delta) then del=delta else del=abs(fid.cosmo.(ind)*delta)
;stop
cosmoh=mk_cosmo(fid,step_param=step_param,step_size=del)
cosmohh=mk_cosmo(fid,step_param=step_param,step_size=(2.d *del))
cosmo_h=mk_cosmo(fid,step_param=step_param,step_size=(-1.d *del))
cosmo_hh=mk_cosmo(fid,step_param=step_param,step_size=(-2.d *del))

clh=mk_cl_tomo(fid,cosmoh,sv)
clhh=mk_cl_tomo(fid,cosmohh,sv)
cl_h=mk_cl_tomo(fid,cosmo_h,sv)
cl_hh=mk_cl_tomo(fid,cosmo_hh,sv)

dcldp=(clh.cl-cl_h.cl-(clhh.cl-cl_hh.cl-(2.d*clh.cl)+(2.d*cl_h.cl))/6.d)/2.d/del
;stop
return,dcldp
end


