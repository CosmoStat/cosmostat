function mk_ddldp_4p,fid,sv,step_param

;Written Tom kitching, August 2008
;Adapted from Adam Amara's code mk_dcldp_4p.pro (2006) and mk_dydp_4p.pro
;This routines calculates the 1st derivative of the m(z) using 4 points.
;INPUT:  fiducial (structure, version 0.11.1 and above)
;        sv       (structure, version 0.11.1 and above)
;        step_param: parameter to vary
;OUTPUT: ddldp, either for z values 
;NOTE: Result agrees with mk_dydp_2p.pro for a given survey to less than
;1e-4 percent

if not tag_exist(fid.cosmo,step_param,index=ind) then begin
   print,'The parameter you are trying to differentiate with respect to is not one the main cosmology parameters.'
   stop
   return,'Error: Parameter not used'
endif
if (fid.calc.verbose ge 1) then print,'varying: '+step_param 
if not tag_check(fid.calc,'delta',val=delta) then delta=0.003d

if (fid.cosmo.(ind) le delta) then del=delta else del=fid.cosmo.(ind)*delta

cosmoh  =mk_cosmo(fid,step_param=step_param,step_size=del,/nopk)
cosmohh =mk_cosmo(fid,step_param=step_param,step_size=(2.d *del),/nopk)
cosmo_h =mk_cosmo(fid,step_param=step_param,step_size=(-1.d *del),/nopk)
cosmo_hh=mk_cosmo(fid,step_param=step_param,step_size=(-2.d *del),/nopk)

mzh   = mk_sne(fid,cosmoh,sv)
mzhh  = mk_sne(fid,cosmohh,sv)
mz_h  = mk_sne(fid,cosmo_h,sv)
mz_hh = mk_sne(fid,cosmo_hh,sv)

ddldp=(mzh.mz_bin-mz_h.mz_bin-(mzhh.mz_bin-mz_hh.mz_bin-$
	(2.d*mzh.mz_bin)+(2.d*mz_h.mz_bin))/6.d)/2.d/del

return,ddldp
end