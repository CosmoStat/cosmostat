pro plt_sv, survey, bin=bin, all_bins=all_bins, overplot=overplot,_extra=extra, unnorm=unnorm

;Written by Anais Rassat, July 2008
; ***************// HELP BEGIN //**************
;PURPOSE: Plots redshift distribution of survey
;INPUT: survey: Survey Structure
;OPTIONAL INPUT: all_bins: select if you want to plot all bins. Do not
;                      use with keyword bin. Not compatible with
;                     'hist' mode survey
;                _extra: any keyword accepted by plot. 
;                bin: number of bin you want. Starts at 1. Not
;                      compatible with 'hist' mode surveys.
;                unnorm: normalized by number of galaxies instead of
;                     int dn/dz dz = 1. Only compatible with 'hist'
;                     mode
; Keywords compatible with 'smail' mode: 
;               all_bins, bin, _extra, overplot
; Keywords compatible with 'hist' mode: 
;               unnorm, _extra
; Note: 'hist' mode not yet compatible with redshift distortions

;EXAMPLES
;---------
; Example 1: Plot a 'smail' type overall distribution and overplot the bins:
; > f = set_fiducial()
; > sv = mk_survey(sv, 'sv1')
; > plt_sv, sv, thick = 4
; > plt_sv, sv, /all_bins, /over
;
; Example 2: Plot a 'hist' type overall distribution, normalized by
; the number of galaxies:
; > f = set_fiducial('wfmos')
; > sv = mk_survey(f, 'sv1')
; > plt_sv, sv, /unnorm
; ***************// HELP END //**************

;If want to plot all bins, first force you to plot bin by bin
if keyword_set(all_bins) then bin = 1
;If only want to plot one bin
if keyword_set(bin) then pz = survey.pz_bin[*,bin-1] else pz = survey.pz_tot
;If don't want a specific yrange then set it to [0,1]
if not keyword_set(yrange) then yrange = [0,1]
xtitle = 'z'
ytitle = 'dN(z)/dz'

if survey.dndztype eq 'smail' then begin
   if not keyword_set(overplot) then plot, survey.z, pz, xtitle=xtitle, ytitle=ytitle, _extra=extra, yrange=yrange else oplot, survey.z, pz, _extra=extra
   
   if keyword_set(all_bins) then begin
      for bin = 1, n_elements(survey.z_min)-1 do begin
         oplot, survey.z, survey.pz_bin[*, bin], _extra=extra
      endfor
   endif
   
endif else if survey.dndztype eq 'hist' then begin
if not keyword_set(unnorm) then yrange = [0.,1.] else yrange=[0.,max(survey.n_g)]
   deltaz = survey.zmed_bin[1]-survey.zmed_bin[0]
   z = [survey.zmed_bin[0]-deltaz, survey.zmed_bin, max(survey.zmed_bin)+deltaz]
   ng = [0.0d, survey.n_g, 0.0d]
   ngnorm = ng & ngnorm[*]=0.d0
   for i = 0, n_elements(ng)-1 do ngnorm[i] = ng[i]/total(ng)

   plot,[0],[0],xran=[min(z), max(z)],yran=yrange, xtitle=xtitle, ytitle=ytitle
   if not keyword_set(unnorm) then oplot, z, ngnorm, psym =10, _extra=extra else oplot, z, ng, psym=10, _extra=extra

endif else print, 'Survey type not supported by plt_sv.pro'

end
