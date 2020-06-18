pro plt_sne, sne, cosmo, errors=errors, overplot=overplot, color=color, _extra=extra

;Written by Thomas Kitching, August 2008.
; ***************// HELP BEGIN //**************
;PURPOSE: Plot supernovae m(k) from sne structure 
;INPUT: sne structure, cosmo structure
;OPTIONAL INPUT: errors: sne covariance for a particluar survey
;                overplot: overplot on previous plot
;                _extra: any keyword accepted by plot or oplot

;EXAMPLE: 
;plt_sne, sne, cosmo
;plt_sne, sne, cosmo, errors=sne_cov, /ylog
; ***************// HELP END //**************

;Set some plotting keywords
if not keyword_set(xlog) then xlog = 0
if not keyword_set(ylog) then ylog = 0
if not keyword_set(xstyle) then xstyle = 0
if not keyword_set(ystyle) then ystyle = 0

xtitle = 'redshift, z'
ytitle = 'm(z)'
if not keyword_set(overplot) then begin
   if not keyword_set(errors) then begin
   	plot, cosmo.evol.z, sne.mz+sne.M0+25., xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog, xstyle=xstyle, ystyle=ystyle, _extra=extra
   endif else begin	
	plot, cosmo.evol.z, sne.mz+sne.M0+25., xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog, xstyle=xstyle, ystyle=ystyle, _extra=extra
      	oploterr,sne.z_bin,sne.mz_bin+sne.M0+25.,sqrt(errors.cov_z)
   endelse
endif else begin
   oplot,cosmo.evol.z, sne.mz, color=color, _extra=extra
endelse

end
