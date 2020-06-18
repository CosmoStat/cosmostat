function xgen,x1,x2,NPOINTS=npoints,LOGPLOT=logplot,DOUBLE_PREC=double_prec

; June 95 - Written by A. Refregier
;
; PURPOSE: generate a vector containing equally spaced numbers
; between x1 and x2. This is typically used to generate an x-axis
; vector in a plot.
; INPUT: x1,x2: interval limits
; OPTIONAL INPUT: npoints: number of points (default=100)
;                 lplot: produce spacing adapted for a log plot
;                 double_prec: produce double precision vector
 
if keyword_set(npoints) then np=npoints else np=100
if keyword_set(double_prec) then yy=dindgen(np) else yy=findgen(np)
if keyword_set(logplot) $
  then xx=10.^(yy*(alog10(x2)-alog10(x1))/float(np)+alog10(x1)) $
  else xx=yy*(float(x2)-float(x1))/float(np)+float(x1)

return,xx
end
