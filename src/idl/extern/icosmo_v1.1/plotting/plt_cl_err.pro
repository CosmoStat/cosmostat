pro plt_cl_err,cl_cov,ind,thick=thick,color=color,fill=fill

;Written by Adam Amara - 4th September 2008
;(modifed from plt_cl_err.pro written by A. Refregier)
; ***************// HELP BEGIN //**************
;PURPOSE:plot the error of the lensing power spectrum 
;INPUTS: cl_cov - lensing power spectrum error from mk_cl_cov.pro
;OPTIONAL INPUTS: thick, color - settings for plotting
;                 fill - fill error boxes
;OUTPUTS: plots
; ***************// HELP END //**************

;; re-arrange so that bin order is [smaller,larger]
;bins=bins(sort(bins))
;; find entry that corresponds to selected bin conbination
;ind1=where(cl_cov.zbins(0,*) eq bins(0))
;ind2=where(cl_cov.zbins(1,ind1) eq bins(1))


; declarations
n_b=n_elements(cl_cov.l_b)
l_l=cl_cov.l_lb & l_h=cl_cov.l_hb & l_m=cl_cov.l_b
cl_m=cl_cov.cl_b(*,ind)
cl_er=cl_cov.cl_err_b(*,ind)
if not keyword_set(color) then color=!p.color

; plot errors
for i=0,n_b-1 do begin
  x=[l_l(i),l_h(i),l_h(i),l_l(i),l_l(i)]>10.^(!x.crange(0))<10.^(!x.crange(1))
  y=([cl_m(i)-cl_er(i),cl_m(i)-cl_er(i),cl_m(i)+cl_er(i),$
          cl_m(i)+cl_er(i),cl_m(i)-cl_er(i)]*$
         l_m(i)*(l_m(i)+1.)/(2.*!pi))>10.^(!y.crange(0))<10.^(!y.crange(1))
  if keyword_set(fill) then polyfill,x,y,color=color $
    else oplot,x,y,thick=thick,color=color
endfor



end

