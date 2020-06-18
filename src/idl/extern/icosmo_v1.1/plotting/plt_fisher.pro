pro plt_fisher,fisher,fisher2,fisher3,fisher4,charsize=charsize,format=format

; Aug 08 - modified header
; Jun 05 - format keyword added by AR
; Mar 03 - Modified by AR to allow plotting of several fisher matrices
; Jan 03 - Modified by AR to plot an arbitrary fisher matrix
; May 01 - Written by A. Refregier
;
; ***************// HELP BEGIN //**************
; PURPOSE: plot the confidence contours for a general fisher matrix
; with an arbitrary number of parameters
; INPUT: fisher: fisher matrix structure given by mk_fisher.pro
; OPTIONAL INPUT: fisher2,3,4: additional fisher matrices which are also
;                     plotted if supplied. This fisher matrix is
;                     assumed to have the same parmeters and central values
;                 format: plot each pair of parameters on a single page 
;			  (default: plot all on the same page) 
; KEYWORDS: charsize: character size for axis labels (default:2)
; OUTPUT: plot
; ***************// HELP END //**************

; declarations
n_p=n_elements(fisher.p)    ; number of parameters
;n_pairs=factorial(n_p)/factorial(2)/factorial(n_p-2) 
;                            ; number of pairs of non-identical parameters
;n_frames=ceil(sqrt(n_pairs)) ; number of frames on a side
if n_p gt 1 and not keyword_set(format) then !p.multi=[0,n_p,n_p]
orig=!p.charsize
if keyword_set(charsize) then !p.charsize=charsize else !p.charsize=2.
col=[4,2,3,5]

; plot contours for every pair of parameters by marginalising over
; the other parameters
pmarg=intarr(n_p)
for i=0,n_p-1 do begin
  for j=0,n_p-1 do begin
    pmarg=replicate(0,n_p)
    pmarg(i)=1 & pmarg(j)=1       ; pick parameters of interest
    margin_fisher,fisher,fisher_m,pmarg
    if n_params() gt 1 then margin_fisher,fisher2,fisher2_m,pmarg
    if n_params() gt 2 then margin_fisher,fisher3,fisher3_m,pmarg
    if n_params() gt 3 then margin_fisher,fisher4,fisher4_m,pmarg
    if i le j then begin
      if i eq j then begin
         plt_fisher_1p,fisher_m,col=col(0)
         if n_params() gt 1 then plt_fisher_1p,fisher2_m,col=col(1),/over
         if n_params() gt 2 then plt_fisher_1p,fisher3_m,col=col(2),/over
         if n_params() gt 3 then plt_fisher_1p,fisher4_m,col=col(3),/over
      endif else begin
         plt_fisher_2p,fisher_m,/central,/nofill,col=col(0)
         if n_params() gt 1 then $
           plt_fisher_2p,fisher2_m,col=col(1),/over,/nofill
         if n_params() gt 2 then $
           plt_fisher_2p,fisher3_m,col=col(2),/over,/nofill
         if n_params() gt 3 then $
           plt_fisher_2p,fisher4_m,col=col(3),/over,/nofill
      endelse
    endif else begin
      plot,[0],[0],xstyle=4,ystyle=4    ; skip this frame
    endelse    
  endfor
endfor
!p.multi=0
!p.charsize=orig 

end
