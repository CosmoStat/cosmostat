pro plt_fisher_1p,fisher,xrange=xrange,over=over,color=color

; Jan 03 - Written by AR
;
; PURPOSE: plot the confidence contours corresponding to Fisher
; matrix for 1 cosmological parameter.
; INPUT: fisher: fisher matrix structure given by mk_fisher.pro
; OPTIONAL INPUT: xran: restricted parameter range for plotting
;                 over: overplot if set
;                 color: color for pdf
; OUTPUT: plot

; declarations
if n_elements(fisher.p) ne 1 then begin
  print,'plt_fisher_1p: fisher matrix must have 1 parameter'
  stop
end
if not keyword_set(color) then color=4
thick=4.
tek_color
sigma=sqrt(fisher.cov(0))
x0=fisher.p(0)
if not keyword_set(xrange) then xrange=4.*sigma*[-1.,1.]+x0

; compute PDF
x=xgen(xrange(0),xrange(1))
p=exp(-(x-x0)^2/(2.*sigma^2))/sqrt(2.*!pi)/sigma

; plot PDF
if not keyword_set(over) then $
  plot,x,p,/nodata,title='!6',$
     xtitle=fisher.pname(0),ytitle='Prob',$
     xrange=xrange,/xstyle
oplot,x,p,thick=thick,color=color

end
