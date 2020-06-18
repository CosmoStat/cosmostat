pro plt_pk, cosmo,z=z, linear=linear, overplot=overplot, color=color, _extra=extra

;Written by Anais Rassat, July 2008.
;July 2008: Modified by An.Rassat: added color keyword.
;PURPOSE: Plot power spectrum P(k) from cosmo structure 
;INPUT: cosmo structure
;OPTIONAL INPUT: z: redshift of power spectrum
;                linear: if set, only linear power spectrum is set
;                overplot: overplot on previous plot
;                _extra: any keyword accepted by plot or oplot

;EXAMPLE: 
;plt_pk, sv
;plt_pk, sv, z=2, linestyle=2, ,/linear,/over

;Setting the redshift
if not keyword_set(z) then pk = create_struct(cosmo.pk) else pk = get_pk(cosmo, z=z)

;Deciding on linear or non-linear
If keyword_set(linear) then pk = pk.pk_l else pk = pk.pk

;Set some plotting keywords
if not keyword_set(xlog) then xlog = 1
if not keyword_set(ylog) then ylog = 1
if not keyword_set(xstyle) then xstyle = 1
if not keyword_set(ystyle) then ystyle = 1
;if not keyword_set(yrange) then yrange = [min(pk), max(pk)*1d1]

xtitle = 'k [h Mpc!U-1!N]'
ytitle = 'P(k) [h!U-3!NMpc!U3!N]'
if not keyword_set(overplot) then begin
   plot, cosmo.pk.k, pk, xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog, xstyle=xstyle, ystyle=ystyle, yrange=yrange, _extra=extra , /nodata
   oplot, cosmo.pk.k, pk, color=color, _extra=extra
endif else begin
   oplot, cosmo.pk.k, pk, color=color, _extra=extra
endelse

end
