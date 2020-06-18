;function mk_dydp_4p,fiducial,cosmo,survey=sv, step_param
function mk_dbaodp_4p,fid,sv,step_param

;Jul 08 - modified by AA made consistent with v0.11
;Written Anais Rassat, July 2008
;Adapted from Adam Amara's code mk_dcldp_4p.pro (2006) and mk_dydp_4p.pro
;This routines calculates the 1st derivative of the tangential BAO scale
;y = r(z)/s using 4 points.
;INPUT:  fiducial (structure, version 0.10.1 and above)
;        cosmo    (structure, version 0.10.1 and above)
;        step_param: parameter to vary
;OPTIONAL INPUT: survey=sv (version 0.10.1 and above)
;OUTPUT: dydp, either for z values given by cosmo.evol.z or if survey
;is specified for sv.z_min

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

;calculate cosmology
cosmoh=mk_cosmo(fid,step_param=step_param,step_size=del,/nopk)
cosmohh=mk_cosmo(fid,step_param=step_param,step_size=2.d*del,/nopk)
cosmo_h=mk_cosmo(fid,step_param=step_param,step_size=(-1.d *del),/nopk)
cosmo_hh=mk_cosmo(fid,step_param=step_param,step_size=(-2.d *del),/nopk)

;calculate the transverse ruler r(z)/s (includes sound horizon s)
baoh=mk_bao(cosmoh,sv)
baohh=mk_bao(cosmohh,sv)
bao_h=mk_bao(cosmo_h,sv)
bao_hh=mk_bao(cosmo_hh,sv)

;calculate derivative
;dydp=(baoh.y-bao_h.y-(baohh.y-bao_hh.y-2.d*bao.h+2.d*y_h)/6.d0)/(2.d*del)

dydp=(baoh.y-bao_h.y-(baohh.y-bao_hh.y-(2.0d*baoh.y)+$
     (2.0d*bao_h.y))/6.0d)/2.0d/del

dydp_bin=(baoh.y_bin-bao_h.y_bin-(baohh.y_bin-bao_hh.y_bin-(2.0d*baoh.y_bin)+$
        (2.0d*bao_h.y_bin))/6.0d)/2.0d/del

dyprimedp=(baoh.yprime-bao_h.yprime-(baohh.yprime-bao_hh.yprime-(2.0d*baoh.yprime)+$
          (2.0d*bao_h.yprime))/6.0d)/2.0d/del

dyprimedp_bin=(baoh.yprime_bin-bao_h.yprime_bin-(baohh.yprime_bin-$
              bao_hh.yprime_bin-(2.0d*baoh.yprime_bin)+(2.0d*bao_h.yprime_bin))/6.0d)/2.0d/del

;construcut a strucutre
dbaodp={dydp:dydp,dyprimedp:dyprimedp,dydp_bin:dydp_bin,dyprimedp_bin:dyprimedp_bin}

return,dbaodp
end


