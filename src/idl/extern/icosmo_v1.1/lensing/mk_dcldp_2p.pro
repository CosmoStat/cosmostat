function mk_dcldp_2p,fid,sv,step_param,keywords=keywords
;Written by Adam Amara 10 May 2006.
;This routines calculates the 1st derivative of a function using 4
;points.
;NEED INSTRUCTIONS !!!!


if not tag_exist(fid,step_param,index=ind) then begin
   print,'The parameter you are trying to differentiate with respect to is not one the main cosmology parameters.'
   return,'Error: Parameter not used'
endif
if (fid.calc.verbose ge 1) then print,'varying: '+step_param 
;if not tag_check(keywords,'delta',val=delta) then delta=0.003d
delta=fid.calc.delta
if (fid.cosmo.(ind) le delta) then del=delta else del=fid.cosmo.(ind)*delta

cosmoh=mk_cosmo(fid,step_param=step_param,step_size=del)
cosmo_h=mk_cosmo(fid,step_param=step_param,step_size=(-1.d *del))

clh=mk_cl_tomo(fid,cosmoh,sv)
cl_h=mk_cl_tomo(fid,cosmo_h,sv)

dcldp=(clh.cl-cl_h.cl)/(2.d*del)
return,dcldp
end


