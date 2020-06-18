pro plt_bao, bao, rad=rad, _extra=extra, overplot=overplot, err=err
; Oct 2008, Cleaned up by An. R.; Removed 'nobin' keyword, added
; plotting of errors.
; Aug 2008, Modified by Anais Rassat. Fixed ytitles.
; Aug 2008, Written by Anais Rassat
; ***************// HELP BEGIN //**************
; PURPOSE: Plot the BAO radial or tangential scales, with or without errors.
; Default is tangential scale, i.e.: r(z)/s or D_A(z)/s
; Also plot radial scale, i.e.: r'(z)/s or c/(H(z)*s)
; INPUT:           bao -   structure output from mk_bao(cosmo, sv)
; OPTIONAL INPUT:  rad - keyword to plot radial modes (default  is
;                        transverse)
;                  err - structure containing bao_cov output from
;                        mk_bao_cov(cosmo, bao, sv)
;                  overplot - if selected, does not create a new plot
;                             but overplots on previous
;                  _extra - any keyword accepted by plot.  If extra
;                           keywords are not accepted by plot, then
;                           these are 'quietly ignored
; OUTPUT: None
; OPTIONAL OUTPUT: None
;
; OBSOLETE KEYWORDS: 
;                  nobin - keyword to plot unbinned modes (binned modes
;                          are survey specific)
;
;----
; Example 1:
;        > fid = set_fiducial()
;        > cosmo = mk_cosmo(fid)
;        > sv = mk_survey(f, 'sv1')
;        > bao = mk_bao(cosmo, sv)
;        > plt_bao, bao
;
; Example 2: 
;        > fid = set_fiducial()
;        > cosmo = mk_cosmo(fid)
;        > sv = mk_survey(f, 'sv1')
;        > bao = mk_bao(cosmo, sv)
;        > bao_cov = mk_bao_cov(cosmo, bao, sv)
;        > plt_bao, bao, err=bao_cov, /rad
; ***************// HELP END //**************
 
tek_color
;tangential and not binned
  x = bao.z
  xtitle = 'z'
  if not keyword_set(rad) then begin ; for transverse scale
     y = bao.y 
     ytitle = 'Transverse BAO Scale: D!DA!N(z)/s'
  endif else begin              ; for radial scale
     y = bao.yprime
     ytitle = 'Radial BAO Scale: c/(H(z)s)'
  endelse
;  endif else begin              ; for binned data
  xtitle = 'z'
  xbin = bao.z_bin
  if not keyword_set(rad) then begin
     ybin = bao.y_bin 
     ytitle = 'Transverse BAO Scale: D!DA!N(z)/s'
; ERRORS NOT YET IMPLEMENTED
        if keyword_set(err) then dy = err.delta_trans
  endif else begin
     ybin = bao.yprime_bin
     ytitle = 'Radial BAO Scale: c/(H(z)s)'
;ERRORS NOT YET IMPLEMENTED
        if keyword_set(err) then dy = err.delta_rad
  endelse
  
;plot or overplot depending on options
  if keyword_set(overplot) then  begin
     oplot, x, y, _extra=extra
  endif else begin
     plot, x, y, xtitle=xtitle, ytitle=ytitle, _extra=extra
     oplot, xbin, ybin, psym = 1,col=2
     
  endelse
  
;overplot errors if desired
  if keyword_set(err) then errplot, xbin, ybin-dy, ybin+dy, col=2
  
end
