pro master_make_pq, delta, lmax, matp, matq, lbins, bintabmax=bintabmax, bintabmin=bintabmin, csbin=csbin
;+
;
; NAME:
;       MASTER_MAKE_PQ
;
; PURPOSE:
;       make binning matix for binning K matrix in master procedure
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       MASTER_MAKE_PQ, delta, lmax, matp, matq
;
; INPUTS:
;       
;       LMAX:   l value at the end of bins
;       DELTA:  width of the bins (first bin starts at l=2, last one 
;                                       ends at lmax)
;
; INPUT KEYWORDS:
;
;       BINTABMAX: array containing ell max of the bins
;       BINTABMIN: array containing ell min of the bins
;    OPTIONNAL:
;       csbin: set to choose a 2l+1 weighting in each bin
;
; OUTPUTS:
;
;       MATP: Pbl matrix
;       MATQ: Qlb matrix
;       LBINS : ell of the centre of each bin
;      
;
; COMMON BLOCKS:
;
;
; SIDE EFFECT:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; MODIFICATION HISTORY:
;
;       Lancien P., 06/07/01
;       Douspis M., 10/06/01
;       Keyword BINTABMAX, BINTABMIN now replace BINTAB
;               CSBIN added, Patanchon G. 03/2004
;       
;       
;-


;nbin            = lmax/delta
;llim            = findgen(nbin+1)*delta
if keyword_set(bintabmin) then begin
 nbin=n_elements(bintabmin)
 llmin=bintabmin
 llmax=bintabmax
end

ell             = lindgen(lmax+1)


MATP  = dblarr(nbin, lmax+1)
MATQ  = dblarr(lmax+1, nbin)

FOR i=0L, nbin-1 DO BEGIN
    FOR ll=llmin(i),llmax(i) DO BEGIN
        IF (KEYWORD_SET(csbin) eq 0) THEN BEGIN
           MATP(i,ll) = 1.d0/2.d0/!dpi*(ll*(ll+1.d0))/double(llmax(i)+1-llmin(i))
           MATQ(ll,i) = 2.*!dpi/double(ll*(ll+1.))
        ENDIF ELSE BEGIN
           MATP(i,ll) = (2.*ll+1.)/total(2.*(findgen(llmax(i)+1-llmin(i))+llmin(i))+1.)
           MATQ(ll,i) = 1.
; print, 'ok'
        ENDELSE
    ENDFOR
ENDFOR

lbins=fltarr(nbin)
for i=0,nbin-1 do lbins(i)=mean([llmin(i),llmax(i)+1])

RETURN

END
