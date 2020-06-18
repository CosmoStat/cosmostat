;+
; NAME:
;        df_del_baseline
;
; PURPOSE:
;   Remove the baseline from a spectrum or an array of spectra
;
; CALLING:
;     TabSpectra_noBaseline =  df_del_baseline(TabSpectra, RMS=RMS,
;     nscale=nscale, baseline = baseline)
;
; INPUTS:
;     TabSpectra -- IDL 1D or 2D array :  TabSpectra[*, 0:N-1]   N input spectra
;    
; OUTPUTS:
;     TabSpectra_noBaseline --   IDL 1D or 2D array :  TabSpectra_noBaseline[*, 0:N-1]   N spectra without baseline
;
; INPUT KEYWORDS:
;  nscaale -- scalar : number of scales used in the decomposition. Default is 6.
;
; EXAMPLE:
;       Remove the baseline of all spectra
;               TabSpectra_noBaseline =  del_baseline(TabSpectra, RMS=RMS) 
;         
; HISTORY:
;	Written: Daniel Machado & Jean-Luc Starck, 2013
;	Sept, 2013 File creation

;--------------------------------------------------------------------------------------------------------
 
function  df_del_baseline, Spectra,  nscale=nscale , rms=rms, baseline=baseline

if N_PARAMS() LT 1  then begin 
   Dec =-1
   print, 'CALLING SEQUENCE: Tabspectra_noBaseline =  df_del_baseline(TabSpectra, RMS=RMS, nscale=nscale, baseline=baseline)'
   goto, DONE
end

input_data = Spectra
vs = size(input_data)
if vs[0] EQ 1 then nspecs=1 $
else nspecs =vs[2]
Nx = vs[1]
errorcurve=0
if keyword_set(RMS) then errorcurve=RMS
if not keyword_set(nsigma) then nsigma = 2. ; detection level in the FDR; nsigma=2 ==> alpha_FDR = 0.05, i.e. 5% of detections may be false detections.

; Option strings for the various denoising steps

if not keyword_set(nscale) then nscale=fix(alog(Nx)) - 1
optscale = '-n ' +  STRC(nscale)
BaselineNscale = nscale + 1

; The first filtering is only used to remove the strongest lines in order to be able to estimate the continuum.
optstring1 = optscale +  ' -f3 -t11 -k -K -P -i20 -s4' ; hard denoising to isolate strong lines

mr1d_filter, input_data, filtered_data, rmsnoise=errorcurve, OPT=optstring1

contnoise = input_data - filtered_data ; creates spectra that are predominantly noise+continuum
continua = 0d * input_data             ; creates a blank matrix the same size as your inputdata.
if nspecs GT 1 then begin
   for p=0, (nspecs-1) do continua[*,p] = (star1d(contnoise[*,p],nscale=BaselineNscale))[*, BaselineNscale-1]
                                ; star1d outputs a matrix with each of
                                ; the scales, per input vector, we
                                ; only want the last  one which is the cleaner continuum estimate.
end else continua = (star1d(contnoise,nscale= BaselineNscale))[*, BaselineNscale-1]

contsubdata = input_data - continua ; new matrix without continua, (but still with noise).
baseline = continua
DONE:

return, contsubdata

END



