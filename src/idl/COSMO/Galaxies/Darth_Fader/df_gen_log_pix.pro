function df_gen_log_pix, lammin, lammax, R=R, delta=delta

if not keyword_set(R) and not keyword_set(delta) then begin
   print, 'No resolution set! Returning...'
   goto, skipall
endif

if keyword_set(R) then begin
   rp = 2.*sqrt(2.*alog(2.))*r
   delta = alog10((rp+0.5)/(rp-0.5)) 
endif else begin
   R=0.5*(10^delta+1)/(10^delta-1)/(2.*sqrt(2.*alog(2.)))
endelse

npix = round(alog10(lammax/lammin)/delta)

logwl = alog10(lammin) + dindgen(npix)*delta

wl = 10.^logwl

skipall:return, wl
end
