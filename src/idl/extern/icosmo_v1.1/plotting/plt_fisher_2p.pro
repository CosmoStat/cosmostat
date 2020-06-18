pro plt_fisher_2p,fisher,xrange=xrange,yrange=yrange,over=over,color=color,$
  prob=prob,nofill=nofill,central=central,_extra=extra

; Oct 08 - Modified by Anais Rassat. _extra keyword passed to plot
;          (before was only passed to oplot).
; Jan 03 - Modified by AR
; May 01 - Written by A. Refregier
;
; PURPOSE: plot the confidence contours corresponding to a 2x2 Fisher
; matrix for 2 cosmological paramaters. This version allows one
; to input contour levels.
; INPUT: fisher: fisher matrix structure given by mk_fisher.pro
; OPTIONAL INPUT: xran,yran: restricted parameter range for plotting
;                 over: overplot if set
;                 col: color for contours (dimension must match number of
;                      contours set by prob; default:4)
;                 prob: probability levels for the contour(s) (default: 68.3%)
;                       (eg. [.683,.955])
;                 nofill: do not fill contours but only draw outlines
;                 central: plot central value as a cross. The value of
;                          the keyword provides the color of the cross.
; OUTPUT: plot

; declarations
if n_elements(fisher.p) ne 2 then begin
  print,'plt_fisher_2p: fisher matrix must have 2 parameters'
  stop
end
if not keyword_set(prob) then prob=.683
n_cont=n_elements(prob)
if not keyword_set(color) then color=replicate(4,n_cont)
thick=4.
tek_color
n_p=1000

; compute contours
v=fisher.f     ; must be 2x2
phi=xgen(0.,2*!pi*1.1,np=n_p)
;chi2=[chisqr_cvf(1.-.683,2),chisqr_cvf(1.-.955,2)]
chi2=fltarr(n_cont)
for i=0,n_cont-1 do chi2(i)=chisqr_cvf(1.-prob(i),2)
r=fltarr(n_cont,n_p)
for i=0,n_cont-1 do begin
  r(i,*)=sqrt(chi2(i)/$
  (v(0,0)*cos(phi)^2+v(1,1)*sin(phi)^2+2.*v(0,1)*cos(phi)*sin(phi)))
endfor

; plot contours
temp=max(prob,cont_max)
if not keyword_set(over) then begin
  plot,r(cont_max,*)*cos(phi)+fisher.p(0),r(cont_max,*)*sin(phi)+fisher.p(1),$
     /ynozero,$
     xtitle=fisher.pname(0),ytitle=fisher.pname(1),$
     xrange=xrange,yrange=yrange,/nodata, _extra=extra
endif
for i=0,n_cont-1 do begin
  x=r(i,*)*cos(phi)+fisher.p(0) & y=r(i,*)*sin(phi)+fisher.p(1)
  if keyword_set(nofill) then oplot,x,y,col=color(i),thick=thick,_extra=extra $
    else polyfill,x,y,noclip=0,color=color(i),_extra=extra
endfor

; place a cross at the central value if requested
if keyword_set(central) then $
  oplot,[fisher.p(0)],[fisher.p(1)],psym=1,col=central

end
