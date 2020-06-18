pro plt_cl, cl,bins,cl_cov,overplot=overplot, color=color,_extra=extra,errors=errors

;Written by Adam Amara - 26th August 2008 (based on the routine
;                        plt_cosmo.pro written by Anais Rassat) 
; ***************// HELP BEGIN //**************
;PURPOSE: Plots any lensing correlation functions
;INPUTS: cl - strucuture output by mk_cl_tomo
;        bins - two element vector of the bins which are
;               correlated. bin index starts with 0 (lowest redshift
;               bin) to nbin-1 (the highest redshift bin), where nbin
;               is the number of bins.
;OPTIONAL INPUTS: overplot - if selected, does not create a new plot
;                            but overplots on previous. 
;                 _extra - any keyword accepted by plot. If extra
;                          keywords are not accepted by plot, then
;                          these are 'quietly ignored'.  
; ----
; Example 1: plt_cl,cl,[0,3]
; This example plots the lensing correlation function between bins 0
; and 3. The plot will use a solid black curve  
;
; Example 2: plt_cl,cl,[1,1],color=2,/overplot,linestyle=2
; This example plots the auto-correlation function of the lensing in
; bin 1. The plots will be produced with red dashed line.  The plot
; will also not erase what is currently displayed in the window.
; i.e. if example 1 is followed by example 2 you will have produced a
; plot with solid black curve showing the [0,3] correlation and a
; dashed green curve with the [1,1] correlation.
;
; ***************// HELP END //**************
on_error,2
tek_color
error_flag=0
; Check that bin ranges are valid
if (n_elements(bins) ne 2) then error_flag=1
if (min(bins) lt 0) then error_flag=1
if (max(bins) gt cl.n_zbin-1) then error_flag=1

; re-arrange so that bin order is [smaller,larger]
bins=bins(sort(bins))
; find entry that corresponds to selected bin conbination
ind1=where(cl.zbins(0,*) eq bins(0))
ind2=ind1(where(cl.zbins(1,ind1) eq bins(1)))
;stop
if not keyword_set(overplot) then begin
   plot,cl.l,cl.l*(cl.l+1.)*cl.cl(*,ind2)/(2.*!pi),/xtype,/ytype,/xstyle,$
        xtitle='!6l',ytitle='l(l+1)C!il!n/(2!7p!6)',title='!6',/nodata,_extra=extra
endif

if keyword_set(errors) then plt_cl_err,cl_cov,ind2,color=2
   oplot,cl.l,cl.l*(cl.l+1.)*cl.cl(*,ind2)/(2.*!pi),color=color,_extra=extra

;if keyword_set(overplot) then begin
;   oplot,cl.l,cl.l*(cl.l+1.)*cl.cl(*,ind2)/(2.*!pi),color=color,_extra=extra
;endif else begin
;   plot,cl.l,cl.l*(cl.l+1.)*cl.cl(*,ind2)/(2.*!pi),/xtype,/ytype,/xstyle,$
;        xtitle='!6l',ytitle='l(l+1)C!il!n/(2!7p!6)',title='!6',/nodata
;   oplot,cl.l,cl.l*(cl.l+1.)*cl.cl(*,ind2)/(2.*!pi),color=color,_extra=extra
;endelse

   print,ind2
End
